import re
import os
import numpy as np
import fitz
from PIL import Image

import re
import os
import numpy as np
import fitz
from PIL import Image

def extract_text_image_data_from_pdf_paths(pdf_paths, extract_images=True):
    """
    Extracts text and image data from PDF files.

    Args:
        pdf_paths (list): A list of strings containing file paths to PDF documents.
        extract_images (bool optional): Whether to extract images from PDF documents.
        
    Returns:
        dict: 

    Raises:
        TypeError: 
        ValueError: 
    """
    
    def extract_text_from_pdf(pdf_path):
        """
        Extracts text from a PDF document, focusing on pages that are most likely to have synthetic content.
        
        The function scans PDF pages for synthetic procedures or other chemical-related information. It uses 
        heuristic-based text processing to clean and identify potentially relevant content. It also identifies 
        consecutive pages that may contain synthetic procedures and performs further text cleaning on those pages.
    
        Args:
            pdf_path (str): Path to the PDF file to be processed.
    
        Returns:
            list: A list of lists, where each sublist contains cleaned synthetic text from consecutive pages 
                  that may contain synthetic procedures.
        """
    
        def find_consecutive_runs_with_padding(nums, num_pages):
            """
            Finds consecutive runs of numbers in a list, adds padding of ±1 to each run, and merges overlapping runs.
    
            Args:
                nums (list): A list of integers representing page numbers. Can be unsorted.
                num_pages (int): Total number of pages in the PDF.
    
            Returns:
                list: A list of lists, where each inner list represents a padded run of consecutive numbers. 
                      Each run is extended by ±1 (if within bounds) and merged with overlapping runs.
            """
            if not nums:
                return []
    
            nums = sorted(nums)
            result = []
            current_run = [nums[0]]
    
            for i in range(1, len(nums)):
                if nums[i] == nums[i - 1] + 1: 
                    current_run.append(nums[i])
                else:
                    result.append(current_run)
                    current_run = [nums[i]]
            result.append(current_run)
    
            padded_result = []
            for run in result:
                start = max(0, run[0] - 1)
                end = min(run[-1] + 1, num_pages - 1)
                padded_result.append(list(range(start, end + 1)))
    
            merged_result = []
            current_run = padded_result[0]
            for i in range(1, len(padded_result)):
                next_run = padded_result[i]
                if current_run[-1] >= next_run[0] - 1:
                    current_run = list(range(current_run[0], next_run[-1] + 1))
                else:
                    merged_result.append(current_run)
                    current_run = next_run
            merged_result.append(current_run)
    
            return merged_result
    
        def clean_chemistry_text(text):
            """
            Cleans text by removing non-relevant lines, including chemical structures, page numbers, 
            and other unwanted content.
    
            Args:
                text (str): The input text from a PDF block.
    
            Returns:
                str: Cleaned text with irrelevant lines removed.
            """
            lines = text.split('\n')
    
            def is_chemical_structure(line):
                """
                Determines if a line represents a chemical structure or formula using heuristics.
    
                Args:
                    line (str): A single line of text.
    
                Returns:
                    bool: True if the line appears to be a chemical structure, False otherwise.
                """
                points = 0
    
                if len(line.strip()) <= 5:
                    points += 1 
    
                chemical_pattern = r'^[A-Z0-9\s\+\-\(\)]+$'
                if re.match(chemical_pattern, line.strip()):
                    points += 1 
    
                chemical_pattern = r"^(?:(?:TBS|TBDMS|TMS|Boc|Cbz|Fmoc|PMB|(?:CF3|F3C)|(?:CCl3|Cl3C)|(?:NO2|O2N)|" \
                                   r"(?:CN|NC)|(?:SO3H|HO3S)|(?:COOH|HOOC)|(?:OH|HO)|(?:NH2|H2N)|Me|Et|Pr|Bu|Ph|Bn|" \
                                   r"Bz|Ac|Ts|Ms|Li|Na|K|Mg|Ca|Fe|Cu|Zn|Pd|Pt|Au|Ag|[HCNOFPSBIK]|Cl|Br|I|R[1-4]|R'?|" \
                                   r"[+-]|[αβ]|(?:cis|trans)|(?:[RS])-|(?:[EZ])-|\d)+)$"
                if re.match(chemical_pattern, line.strip()):
                    points += 3
    
                chemical_chars = set('OHNClBrFPSMet')
                char_count = sum(1 for c in line if c in chemical_chars)
                if len(line.strip()) > 0:
                    chemical_ratio = char_count / len(line.strip())
                    if chemical_ratio > 0.8:
                        points += 1 
    
                return points >= 3
    
            def is_page_number(line):
                """
                Determines if a line is a page number.
    
                Args:
                    line (str): A single line of text.
    
                Returns:
                    bool: True if the line is a page number, False otherwise.
                """
                return re.match(r'^\d+$', line.strip()) is not None
    
            cleaned_lines = []
            for i, line in enumerate(lines):
                if not line.strip():
                    continue
                if is_page_number(line) and (i == 0 or i == len(lines) - 1):
                    continue
                if is_chemical_structure(line):
                    continue
                cleaned_lines.append(line)
    
            return '\n'.join(cleaned_lines)

    
        WORDS_TO_CHECK = [' mg', ' ml', ' mL', ' g', ' mol', ' mmol']
        MUST_INCLUDE = ['yield', 'yielded', 'yields', 'obtain', 'obtained', 'obtains',
                        'synthesis', 'synthesized', 'preparation', 'prepared', 'nmr',
                        'm.p.', 'i.r.', 'm/z', 'purified', 'extracted', 'procedure',
                        'isolated', 'title compound']
        WORDS_TO_CHECK_PATTERN = r'(?i)(?<=\d)(?:' + '|'.join(re.escape(word) for word in WORDS_TO_CHECK) + r')(?=\W|$)'

        if not isinstance(pdf_path, str):
            raise TypeError("The 'pdf_path' argument must be a string.")
        
        if not os.path.exists(pdf_path):
            raise TypeError("The 'pdf_path' argument must be a valid path.")
        
        if not os.path.isfile(pdf_path):
            raise TypeError("The 'pdf_path' argument must be a file.")

        try:
            doc = fitz.open(pdf_path)
        except:
            raise TypeError("The 'pdf_path' argument must be a pdf.")
            
        pages_that_may_contain_synthetic_procedures = []

        all_blocks = []
        for page_num in range(len(doc)):
            page = doc[page_num]
            blocks = page.get_text("blocks")

            elements = []
            for block in blocks:
                if block[6] == 0:
                    elements.append({
                        'type': 'text',
                        'coords': block[0:4],
                        'content': block[4]
                    })

            elements.sort(key=lambda x: x['coords'][1])
            all_blocks.append(elements)
            
            page_content = []
            for element in elements:
                page_content.append(element['content'])
            joined_blocks = '\n'.join(page_content)
    
            if re.search(WORDS_TO_CHECK_PATTERN, joined_blocks.lower()):
                if any(word in joined_blocks.lower() for word in MUST_INCLUDE):
                    pages_that_may_contain_synthetic_procedures.append(page_num)
    
        consecutive_runs = find_consecutive_runs_with_padding(pages_that_mayc_contain_synthetic_procedures, len(doc))
    
        pages_with_synthetic_text = {}
        for consecutive_run in consecutive_runs:
            for page_num in consecutive_run:
                pages_with_synthetic_text[page_num] = all_blocks[page_num]
    
        doc.close()
        return pages_with_synthetic_text

    extracted_pdf_dict = {}
    for i, pdf_path in enumerate(pdf_paths):
        extracted_pdf_dict[i] = extract_text_from_pdf(pdf_path)
    
    return extracted_pdf_dict
