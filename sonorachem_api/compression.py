import io
import zipfile
import tempfile
import base64
from pathlib import Path
from ast import literal_eval

def process_zip_file(zip_path):
    """
    Process a ZIP file and extract JSON data from all .txt files within it.

    Args:
        zip_path (str): The path to the ZIP file.

    Returns:
        dict: The combined JSON data from all processed text files.
    """
    combined_json_data = {}
    with zipfile.ZipFile(zip_path, 'r') as archive:
        for file_name in archive.namelist():
            if file_name.endswith('.txt'):
                text = process_text_file(archive, file_name)
                json_data = literal_eval(text)
                combined_json_data[file_name] = json_data
    return combined_json_data
    
def process_text_file(archive, file_name):
    """
    Process a single text file within a ZIP archive.

    Args:
        archive (zipfile.ZipFile): The ZIP archive object.
        file_name (str): The name of the text file to process.

    Returns:
        str: The concatenated content of the text file.
    """
    with archive.open(file_name) as text_file:
        content = text_file.read().decode('utf-8')
    return ''.join(content.splitlines())

def compress_data(input_data):
    """
    Compresses input data using ZIP compression and encodes it in base64.

    Args:
        input_data: The data to be compressed (can be any type that can be converted to string)

    Returns:
        str: Base64 encoded string of the compressed data
    """
    try:
        memory_file = io.BytesIO()
        
        with zipfile.ZipFile(memory_file, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
            zip_file.writestr('input_data.txt', str(input_data))
        
        memory_file.seek(0)
        zip_content = memory_file.getvalue()
        compressed_data = base64.b64encode(zip_content).decode('utf-8')
        
        return compressed_data
    
    except Exception as e:
        raise Exception(f"Error compressing data: {str(e)}")
    
    finally:
        memory_file.close()

def decompress_data(encoded_zip: str) -> List[Any]:
    """
    Decompresses a base64 encoded ZIP string back into its original data.

    Args:
        encoded_zip (str): Base64 encoded string containing compressed ZIP data

    Returns:
        List[Any]: List containing the decompressed data

    Raises:
        Exception: If there's an error during decompression or file processing
    """
    temp_zip_path = None
    try:
        zip_content = base64.b64decode(encoded_zip)
        
        with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as temp_zip:
            temp_zip.write(zip_content)
            temp_zip_path = temp_zip.name
        
        json_data = _process_zip_file(temp_zip_path)
        output_data = list(json_data.values())
        
        return output_data

    except Exception as e:
        raise Exception(f"Error decompressing data: {str(e)}")
    
    finally:
        if temp_zip_path and Path(temp_zip_path).exists():
            Path(temp_zip_path).unlink()
