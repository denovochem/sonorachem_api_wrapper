import re
import os
import numpy as np
import fitz
from PIL import Image

def extract_text_image_data_from_pdf_paths(pdf_paths):
    """
    Extracts text and image data from PDF files.

    Args:
        pdf_paths (list): A list of strings containing file paths to PDF documents.

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

        def rect_distance(rect1, rect2):
            """
            Calculate the minimum distance between two rectangles.

            Args:
                rect1 (tuple): First rectangle coordinates (x1, y1, x2, y2).
                rect2 (tuple): Second rectangle coordinates (x3, y3, x4, y4).

            Returns:
                float: Minimum distance between rectangles, 0 if overlapping.
            """
            
            x1, y1, x2, y2 = rect1
            x3, y3, x4, y4 = rect2

            # Check for overlap
            if x1 <= x4 and x2 >= x3 and y1 <= y4 and y2 >= y3:
                return 0  # Rectangles overlap

            # Calculate horizontal and vertical distances
            dx = max(x1 - x4, x3 - x2, 0)
            dy = max(y1 - y4, y3 - y2, 0)

            # Euclidean distance
            return (dx**2 + dy**2) ** 0.5

        def should_merge(rect1, rect2, distance_threshold):
            """
            Check if two bounding boxes should be merged based on distance.

            Args:
                rect1 (Rect): First bounding box.
                rect2 (Rect): Second bounding box.
                distance_threshold (float): Maximum distance for merging.

            Returns:
                bool: True if boxes should be merged, False otherwise.
            """
            
            return rect_distance(rect1, rect2) <= distance_threshold

        def merge_rects(rect1, rect2):
            """
            Merge two bounding boxes into one.

            Args:
                rect1 (Rect): First bounding box.
                rect2 (Rect): Second bounding box.

            Returns:
                Rect: New merged bounding box.
            """
            
            x0 = min(rect1.x0, rect2.x0)
            y0 = min(rect1.y0, rect2.y0)
            x1 = max(rect1.x1, rect2.x1)
            y1 = max(rect1.y1, rect2.y1)
            return fitz.Rect(x0, y0, x1, y1)

        def merge_bounding_boxes(bounding_boxes, distance_threshold):
            """
            Iteratively merge bounding boxes until no more merges are possible.

            Args:
                bounding_boxes (list): List of bounding boxes.
                distance_threshold (float): Maximum distance for merging.

            Returns:
                list: List of merged bounding boxes.
            """
            
            merged = True
            while merged:
                merged = False
                new_boxes = []
                used_indices = set()
                for i in range(len(bounding_boxes)):
                    if i in used_indices:
                        continue
                    rect1 = bounding_boxes[i]
                    new_rect = rect1
                    for j in range(i + 1, len(bounding_boxes)):
                        if j in used_indices:
                            continue
                        rect2 = bounding_boxes[j]
                        if should_merge(new_rect, rect2, distance_threshold):
                            new_rect = merge_rects(new_rect, rect2)
                            used_indices.add(j)
                            merged = True
                    new_boxes.append(new_rect)
                    used_indices.add(i)
                bounding_boxes = new_boxes
            return bounding_boxes

        def group_drawings_by_proximity(drawings, distance_threshold=3):
            """
            Group drawings based on proximity.

            Args:
                drawings (list): List of drawing dictionaries.
                distance_threshold (float, optional): Maximum distance for grouping. Defaults to 3.

            Returns:
                list: List of drawing groups.
            """
            
            groups = []
            for drawing in drawings:
                rect = drawing["rect"]
                added_to_group = False
                for group in groups:
                    for other_drawing in group:
                        other_rect = other_drawing["rect"]
                        if rect_distance(rect, other_rect) <= distance_threshold:
                            group.append(drawing)
                            added_to_group = True
                            break
                    if added_to_group:
                        break
                if not added_to_group:
                    groups.append([drawing])
            return groups

        def expand_bbox_with_text(bbox, text_elements, text_distance_threshold=10):
            """
            Expand bounding box to include nearby text.

            Args:
                bbox (tuple): Original bounding box coordinates.
                text_elements (list): List of text elements.
                text_distance_threshold (float, optional): Maximum distance to include text. Defaults to 10.

            Returns:
                Rect: Expanded bounding box.
            """
            
            x0, y0, x1, y1 = bbox
            text_elements_to_merge = []
            for text in text_elements:
                for line in text.get('lines', []):
                    for span in line.get('spans', []):
                        span_bbox = fitz.Rect(span["bbox"])
                        if rect_distance(bbox, span_bbox) <= text_distance_threshold:
                            text_elements_to_merge.append(text)
                            x0 = min(x0, span_bbox.x0)
                            y0 = min(y0, span_bbox.y0)
                            x1 = max(x1, span_bbox.x1)
                            y1 = max(y1, span_bbox.y1)

            return fitz.Rect(x0, y0, x1, y1)

        def get_drawing_bboxes(page, page_num, image_index, return_images=True, bbox_padding=5):
            """
            Generate bounding boxes and images for drawings in the PDF.

            Args:
                page (Page): PDF page object.
                page_num (int): Page number.
                image_index (int): image index number.
                return_images (bool, optional): Whether to return images. Defaults to True.
                bbox_padding (int, optional): Padding for bounding boxes. Defaults to 5.

            Returns:
                list: List of dictionaries containing page number, bbox, and image data.
            """

            returned_data = []
            extracted_images = page.get_image_info()
            for extracted_image in extracted_images:

                returned_image_dict = {
                    'page_num': page_num,
                    'image_index': image_index,
                    'bbox': list(extracted_image['bbox']),
                    'image': np.array([])
                }

                if return_images:
                    pix = page.get_pixmap(clip=list(extracted_image['bbox']))
                    pil_image = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
                    returned_image_dict['image'] = np.array(pil_image)

                returned_data.append(returned_image_dict)
                
                image_index += 1

            drawings = page.get_drawings()

            groups = group_drawings_by_proximity(drawings, distance_threshold=50)

            bounding_boxes = []
            for group in groups:
                x0 = min(d["rect"].x0 for d in group)
                y0 = min(d["rect"].y0 for d in group)
                x1 = max(d["rect"].x1 for d in group)
                y1 = max(d["rect"].y1 for d in group)
                bounding_boxes.append(fitz.Rect(x0, y0, x1, y1))

            merged_bounding_boxes = merge_bounding_boxes(bounding_boxes, distance_threshold=50)

            text_elements = page.get_text("dict")["blocks"]

            expanded_bounding_boxes = [expand_bbox_with_text(bbox, text_elements) for bbox in merged_bounding_boxes]

            for expanded_bounding_box in expanded_bounding_boxes:
                expanded_bounding_box.x0 -= bbox_padding
                expanded_bounding_box.y0 -= bbox_padding
                expanded_bounding_box.x1 += bbox_padding
                expanded_bounding_box.y1 += bbox_padding

                returned_image_dict = {
                    'page_num': page_num,
                    'image_index': image_index,
                    'bbox': [expanded_bounding_box.x0, expanded_bounding_box.y0, expanded_bounding_box.x1, expanded_bounding_box.y1],
                    'image': np.array([])
                }

                if return_images:
                    pix = page.get_pixmap(clip=expanded_bounding_box)
                    pil_image = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
                    returned_image_dict['image'] = np.array(pil_image)

                returned_data.append(returned_image_dict)
                
                image_index += 1

            returned_data = filter_common_bboxes(returned_data, x=3, y=int(np.ceil(len(doc)/2)))
            returned_data = filter_bboxes_aspect_ratio(returned_data, max_ratio=5) 
            returned_data = filter_bboxes_size(returned_data, (page.rect.width, page.rect.height), max_percentage=0.8) 

            return returned_data, image_index

        def is_similar_bbox(bbox1, bbox2, x):
            """
            Check if two bounding boxes are similar within x%.

            Args:
                bbox1 (tuple): First bounding box coordinates.
                bbox2 (tuple): Second bounding box coordinates.
                x (float): Similarity threshold percentage.

            Returns:
                bool: True if boxes are similar, False otherwise.
            """

            x0_diff = abs(bbox1[0] - bbox2[0]) / ((bbox1[0] + bbox2[0]) / 2) * 100
            y0_diff = abs(bbox1[1] - bbox2[1]) / ((bbox1[1] + bbox2[1]) / 2) * 100
            x1_diff = abs(bbox1[2] - bbox2[2]) / ((bbox1[2] + bbox2[2]) / 2) * 100
            y1_diff = abs(bbox1[3] - bbox2[3]) / ((bbox1[3] + bbox2[3]) / 2) * 100

            return x0_diff <= x and y0_diff <= x and x1_diff <= x and y1_diff <= x

        def filter_bboxes_aspect_ratio(returned_data, max_ratio=3):
            """
            Filter out bounding boxes with extreme aspect ratios.
        
            Args:
                returned_data (list): List of bounding box data.
                max_ratio (float): Maximum allowed aspect ratio (width:height or height:width).
        
            Returns:
                list: Filtered list of bounding box data.
            """
            filtered_data = []
            
            for item in returned_data:
                bbox = item['bbox']
                width = bbox[2] - bbox[0]
                height = bbox[3] - bbox[1]
                
                aspect_ratio = width / height if height != 0 else float('inf')
                inverse_aspect_ratio = height / width if width != 0 else float('inf')
                
                if aspect_ratio <= max_ratio and inverse_aspect_ratio <= max_ratio:
                    filtered_data.append(item)
                    
            return filtered_data
        
        def filter_bboxes_size(returned_data, page_size, max_percentage=0.5):
            """
            Filter out bounding boxes that are too large relative to the page size.
        
            Args:
                returned_data (list): List of bounding box data.
                page_size (tuple): Tuple of (page_width, page_height).
                max_percentage (float): Maximum allowed size as a percentage of page area (0.0 to 1.0).
        
            Returns:
                list: Filtered list of bounding box data.
            """
            filtered_data = []
            page_area = page_size[0] * page_size[1]
            
            for item in returned_data:
                bbox = item['bbox']
                bbox_width = bbox[2] - bbox[0]
                bbox_height = bbox[3] - bbox[1]
                bbox_area = bbox_width * bbox_height
                
                area_percentage = bbox_area / page_area
                
                if area_percentage <= max_percentage:
                    filtered_data.append(item)
                    
            return filtered_data

        def filter_common_bboxes(returned_data, x, y):
            """
            Filter out bounding boxes that appear on more than y pages within x% similarity.

            Args:
                returned_data (list): List of bounding box data.
                x (float): Similarity threshold percentage.
                y (int): Maximum number of occurrences.

            Returns:
                list: Filtered list of bounding box data.
            """
            bbox_occurrences = {}

            for item in returned_data:
                bbox = tuple(item['bbox'])
                is_common = False

                for key in bbox_occurrences:
                    if is_similar_bbox(bbox, key, x):
                        bbox_occurrences[key] += 1
                        is_common = True
                        break

                if not is_common:
                    bbox_occurrences[bbox] = 1

            filtered_data = []
            for item in returned_data:
                bbox = tuple(item['bbox'])
                is_common = False

                for key in bbox_occurrences:
                    if is_similar_bbox(bbox, key, x) and bbox_occurrences[key] > y:
                        is_common = True
                        break

                if not is_common:
                    filtered_data.append(item)

            return filtered_data
    
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
        all_images = []
        image_index = 0
        for page_num in range(len(doc)):
            page = doc[page_num]
            
            image_bboxes, image_index = get_drawing_bboxes(page, page_num, image_index)
            all_images.extend(image_bboxes)
            blocks = page.get_text("blocks")

            elements = []
            for block in blocks:
                if block[6] == 0:
                    elements.append({
                        'type': 'text',
                        'coords': block[0:4],
                        'content': block[4]
                    })

            for img in image_bboxes:
                elements.append({
                    'type': 'image',
                    'coords': img['bbox'],
                    'content': f"[###IMAGE_{img['image_index']}###]"
                })

            elements.sort(key=lambda x: x['coords'][1])

            page_content = []
            for element in elements:
                page_content.append(element['content'])

            joined_blocks = '\n'.join(page_content)
            all_blocks.append(joined_blocks)
    
            if re.search(WORDS_TO_CHECK_PATTERN, joined_blocks.lower()):
                if any(word in joined_blocks.lower() for word in MUST_INCLUDE):
                    pages_that_may_contain_synthetic_procedures.append(page_num)
    
        consecutive_runs = find_consecutive_runs_with_padding(pages_that_may_contain_synthetic_procedures, len(doc))
    
        pages_with_synthetic_text = {}
        for consecutive_run in consecutive_runs:
            for page_num in consecutive_run:
                synthetic_text = clean_chemistry_text(all_blocks[page_num])
                pages_with_synthetic_text[page_num] = synthetic_text
    
        doc.close()
        return pages_with_synthetic_text, all_images

    extracted_pdf_dict = {}
    for i, pdf_path in enumerate(pdf_paths):
        extracted_pdf_dict[i] = extract_text_from_pdf(pdf_path)
    
    return extracted_pdf_dict
