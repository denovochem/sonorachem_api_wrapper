import requests
import json
import time
import io
import zipfile
import tempfile
import base64
from pathlib import Path
from ast import literal_eval
from typing import Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*')  
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

from .utils import randomize_smiles, randomize_reaction_smiles

class SonoraChemAPIWrapper:
    """
    Python wrapper for SonoraChem API.
    """

    def __init__(
        self,
        api_key: str,
        base_url: Optional[str] = 'https://hwsvflxtqh.execute-api.us-east-2.amazonaws.com/denovochem_api_stage/run',
    ):
        """
        SonoraChemAPIWrapper constructor.

        Args:
            api_key (str): an API key to access the service.
            base_url (str, optional): base url for the service. If not provided it will default to
                the De Novo Chem AWS server
        """
        self._api_key = api_key
        self._base_url = base_url
        self._headers = self._construct_headers()
        self._validate_api_key()

    def _validate_api_key(self):
        """
        Validate the API key by sending a POST request to the usage URL.
    
        This method attempts to validate the API key stored in the instance variable
        self._api_key. It sends a POST request to the usage URL and checks the response.
    
        Raises:
            ValueError: If the API key validation fails or if there's an unexpected
                        response from the server.
    
        Returns:
            None
        """
        post_request_data = {
            "endpoint": 'get_usage',
        }

        post_request_data = {"input": post_request_data}

        response = self._send_post_request(self._base_url, self._headers, post_request_data)
                
    def _construct_headers(self) -> dict:
        """
        Construct header, required for all requests.

        Returns:
            dict: dictionary containing the "Content-Type" and the
                "x-api-key".
        """
        return {"Content-Type": "application/json", "x-api-key": self._api_key}

    def _send_post_request(self, 
                          url: str, 
                          headers: Dict[str, str], 
                          data: Dict[str, Any], 
                          timeout: int = 120, 
                          verify_ssl: bool = True) -> Optional[Dict[str, Any]]:
        """
        Sends a POST request to the specified URL with given headers and data.
    
        This function attempts to send a POST request and handles various potential
        errors that might occur during the process. It also allows for customization
        of the request timeout and SSL verification.
    
        Args:
            url (str): The URL to send the POST request to.
            headers (Dict[str, str]): A dictionary of HTTP headers to send with the request.
            data (Dict[str, Any]): A dictionary of data to send in the body of the POST request.
            timeout (int, optional): The maximum number of seconds to wait for a response. 
                Defaults to 30 seconds.
            verify_ssl (bool, optional): Whether to verify the server's TLS certificate. 
                Defaults to True.
    
        Returns:
            Optional[Dict[str, Any]]: The JSON-decoded response content if the request was 
            successful and returned valid JSON. Returns None if the request failed or 
            the response wasn't valid JSON.
    
        Raises:
            requests.RequestException: For any request-related errors.
            json.JSONDecodeError: If the response content is not valid JSON.
            ValueError: If the URL is not valid.
            
        """
        try:
            if not url.startswith(('http://', 'https://')):
                raise ValueError("Invalid URL. Must start with 'http://' or 'https://'")
    
            response = requests.post(url, 
                                     headers=headers, 
                                     json=data,  # Use json parameter for automatic JSON encoding
                                     timeout=timeout, 
                                     verify=verify_ssl)
            
            response.raise_for_status()
    
            return response.json()
    
        except requests.RequestException as e:
            print(f"Request failed: {e}")
            if hasattr(e, 'response'):
                print(f"Response status code: {e.response.status_code}")
                print(f"Response content: {e.response.text}")
        except json.JSONDecodeError:
            print("Response was not valid JSON")
        except ValueError as e:
            print(f"Value error: {e}")
        
        return None

    def set_base_url(self, base_url: str) -> None:
        """
        Set base url.

        Args:
            base_url (str): base url for the service to set.
        """
        self._base_url = base_url

    def set_api_key(self, api_key: str):
        """
        Set the API key.

        This method also rebuilds the headers.

        Args:
            api_key (str): an API key to access the service.

        Examples:
            Set an API key:

            >>> sonorachem_api_wrapper.set_api_key('API_KEY')
        """
        logger.info("Set API key to {}".format(api_key))
        self._api_key = api_key
        self._headers = self._construct_headers()
        self._validate_api_key()
        
    def get_api_usage(self):
        """
        Retrieves the API usage statistics for the current API key.
    
        This function sends a POST request to the API's usage endpoint to fetch the
        usage statistics associated with the provided API key. The response is then
        returned as a dictionary.
    
        Returns:
            dict: A dictionary containing the API usage data.
    
        Raises:
            HTTPError: If the request to the API fails.
            ValueError: If the response from the API is not in the expected format.
        """
        post_request_data = {
            "endpoint": 'get_usage',
        }

        post_request_data = {"input": post_request_data}
    
        output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
    
        return output_data

    def is_valid_smiles(self, smiles):
        """
        Sanitizes and checks if the input is a valid SMILES string, potentially with one or more fragments.
        Removes atom mapping and isotopes if present and issues a warning.
    
        Parameters:
        smiles (str): The SMILES string to be checked.
    
        Returns:
        bool: True if the SMILES string is valid, False otherwise.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                return False
            
            if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                print("Warning: Atom mapping found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
            
            if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                print("Warning: Isotopes found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetIsotope(0)
            
            Chem.SanitizeMol(mol)
            
            return True
        except Exception as e:
            print(f"Error: {e}")
            return False

    def is_valid_reaction_smiles(self, reaction_smiles):
        """
        Sanitizes and checks if the input is a valid reaction SMILES.
        Removes atom mapping and isotopes if present and issues a warning.
    
        Parameters:
        reaction_smiles (str): The reaction SMILES string to be checked.
    
        Returns:
        bool: True if the reaction SMILES is valid, False otherwise.
        """
        try:
            reaction = AllChem.ReactionFromSmarts(reaction_smiles)
            
            if reaction is None:
                return False
            
            for mol in reaction.GetReactants():
                if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Atom mapping found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
    
            for mol in reaction.GetProducts():
                if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Atom mapping found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
            
            for mol in reaction.GetReactants():
                if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Isotopes found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetIsotope(0)
    
            for mol in reaction.GetProducts():
                if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Isotopes found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetIsotope(0)
            
            Chem.rdChemReactions.SanitizeRxn(reaction)
            
            if reaction.Validate()[1] != 0:
                return False
            
            return True
        except Exception as e:
            print(f"Error: {e}")
            return False

    def _predict(self, endpoint, input_data, input_data_type='smiles', model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Parent function to make a prediction.

        Args:
            endpoint (str): The endpoint to use for the prediction.
            input_data (str): input data.
            input_data_type (str, optional): The input data type used for input validation. 
                Must be one of 'smiles' or 'rxn_smiles'. Defaults to 'smiles'.
            model_version (str, optional): The version of the model to use. Defaults to 'latest'.
            sampling_method (str, optional): The method used for sampling predictions. 
                Must be one of 'top_k', 'greedy', or 'sampling'. Defaults to 'greedy'.
            seq_length (int, optional): The maximum sequence length for the model input. 
                Defaults to 512.
            beam_size (int, optional): The beam size for beam search (if applicable). 
                Defaults to 5.
            temperature (float, optional): The temperature parameter for controlling randomness 
                in sampling. Must be a positive float. Higher values increase randomness. 
                Defaults to 1.0.

        Returns:
            dict: A dictionary containing the input, output, status, and execution time.

        Raises:
            ValueError: 
              - If `input_data` is not a string.
              - If `model_version` is not a string.
              - If `input_data_type` is not one of 'smiles', 'rxn_smiles'.
              - If `sampling_method` is not one of 'top_k', 'greedy', or 'sampling'.
              - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
              - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
              - If `temperature` is not a positive float.
              - if `input_data` is not a valid SMILES string.
        """
        if not isinstance(input_data, str):
            raise TypeError("The 'input_data' argument must be a string.")
            
        if input_data_type not in ['smiles', 'rxn_smiles']:
            raise ValueError("Invalid 'input_data_type'. Must be 'smiles', 'rxn_smiles'.")

        if not isinstance(model_version, str):
            raise TypeError("The 'model_version' argument must be a string.")
        
        if sampling_method not in ['top_k', 'greedy', 'sampling']:
            raise ValueError("Invalid sampling method. Must be 'top_k', 'greedy', or 'sampling'.")

        if not isinstance(seq_length, int):
            raise TypeError("The 'seq_length' argument must be an integer.")
        if seq_length <= 0 or seq_length > 512:
            raise ValueError("The 'seq_length' argument must be greater than 0 and less than or equal to 512.")
    
        if not isinstance(beam_size, int):
            raise TypeError("The 'beam_size' argument must be an integer.")
        if beam_size <= 0 or beam_size > 16:
            raise ValueError("The 'beam_size' argument must be greater than 0 and less than or equal to 16.")
        
        if beam_size == 1:
            sampling_method = 'greedy'
        
        if temperature <= 0:
            raise ValueError("The 'temperature' argument must be a positive float.")

        if input_data_type == 'smiles':
            valid_smiles = self.is_valid_smiles(input_data)
            if not valid_smiles:
                raise ValueError("The 'input_data' argument is not a valid SMILES string.")
                
        if input_data_type == 'rxn_smiles':
            valid_rxn_smiles = self.is_valid_reaction_smiles(input_data)
            if not valid_rxn_smiles:
                raise ValueError("The 'input_data' argument is not a valid reaction SMILES string.")

        post_request_data = {
            "endpoint": endpoint,
            "data": {
                "model_version": model_version,
                "input_data": input_data,
                "kwargs": {
                    "sampling_method": sampling_method,
                    "seq_length": seq_length,
                    "beam_size": beam_size,
                    "temperature": temperature
                }
            }
        }

        post_request_data = {"input": post_request_data}

        start = time.time()
        output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
        returned_data = {
            'input': post_request_data['input'],
            'output': output_data['output'],
            'status': output_data['status'],
            'execution_time': time.time() - start
        }
    
        return returned_data

    def _batch_predict(self, endpoint, input_data,  input_data_type='smiles', model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Parent function to make batch predictions.

        Args:
            endpoint (str): The endpoint to use for the prediction.
            input_data (list of str): A list of SMILES strings.
            input_data_type (str, optional): The input data type used for input validation. 
                Must be one of 'smiles' or 'rxn_smiles'. Defaults to 'smiles'.
            model_version (str, optional): The version of the model to use. Defaults to 'latest'.
            sampling_method (str, optional): The method used for sampling predictions. 
                Must be one of 'greedy', or 'sampling'. Defaults to 'greedy'.
            seq_length (int, optional): The maximum sequence length for the model input. 
                Defaults to 512.
            beam_size (int, optional): The beam size for beam search (if applicable). 
                Defaults to 5.
            temperature (float, optional): The temperature parameter for controlling randomness 
                in sampling. Must be a positive float. Higher values increase randomness. 
                Defaults to 1.0.
            batch_size (int, optional): The batch size for batch inference. Defaults to 64.

        Returns:
            dict: A dictionary containing the input, output, status, and execution time.

        Raises:
            ValueError: 
              - If `input_data` is not a list of strings.
              - If `input_data_type` is not one of 'smiles', 'rxn_smiles'.
              - If `model_version` is not a string.
              - If `sampling_method` is not one of 'greedy', or 'sampling'.
              - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
              - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
              - If `temperature` is not a positive float.
              - if `batch_size` is not an integer, or if it is not greater than 0 or exceeds 256.
              - If any element in `input_data` is not a valid SMILES string.
        """
        if not isinstance(input_data, list) or not all(isinstance(item, str) for item in input_data):
            raise TypeError("The 'input_data' argument must be a list of strings.")

        if input_data_type not in ['smiles', 'rxn_smiles']:
            raise ValueError("Invalid 'input_data_type'. Must be 'smiles', 'rxn_smiles'.")

        if not isinstance(model_version, str):
            raise TypeError("The 'model_version' argument must be a string.")

        if sampling_method not in ['greedy', 'sampling']:
            raise ValueError("Invalid sampling method. Must be 'greedy', or 'sampling'.")

        if not isinstance(seq_length, int):
            raise TypeError("The 'seq_length' argument must be an integer.")
        if seq_length <= 0 or seq_length > 512:
            raise ValueError("The 'seq_length' argument must be greater than 0 and less than or equal to 512.")

        if not isinstance(beam_size, int):
            raise TypeError("The 'beam_size' argument must be an integer.")
        if beam_size <= 0 or beam_size > 16:
            raise ValueError("The 'beam_size' argument must be greater than 0 and less than or equal to 16.")

        if beam_size == 1:
            sampling_method = 'greedy'

        if temperature <= 0:
            raise ValueError("The 'temperature' argument must be a positive float.")

        if not isinstance(batch_size, int):
            raise TypeError("The 'batch_size' argument must be an integer.")
        if batch_size <= 0 or batch_size > 256:
            raise ValueError("The 'batch_size' argument must be greater than 0 and less than or equal to 256.")

        if input_data_type == 'smiles':
            for smiles in input_data:
                valid_smiles = self.is_valid_smiles(smiles)
                if not valid_smiles:
                    raise ValueError(f"The SMILES string '{smiles}' is not valid.")
                
        if input_data_type == 'rxn_smiles':
            for rxn_smiles in input_data:
                valid_rxn_smiles = self.is_valid_reaction_smiles(rxn_smiles)
                if not valid_rxn_smiles:
                    raise ValueError(f"The reaction SMILES string '{rxn_smiles}' is not valid.")

        post_request_data = {
            "endpoint": endpoint,
            "data": {
                "model_version": model_version,
                "input_data": input_data,
                "kwargs": {
                    "sampling_method": sampling_method,
                    "seq_length": seq_length,
                    "beam_size": beam_size,
                    "temperature": temperature,
                    "batch_size": batch_size
                }
            }
        }

        post_request_data = {"input": post_request_data}

        start = time.time()
        output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
        returned_data = {
            'input': post_request_data['input'],
            'output': output_data['output'],
            'status': output_data['status'],
            'execution_time': time.time() - start
        }
    
        return returned_data

    def predict_procedures_retro_template_free(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict retrosynthetic procedures for a given SMILES string using a template-free approach.
        """
        return self._predict("procedures_retro_template_free", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_procedures_retro_template_free(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict retrosynthetic procedures for SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_procedures_retro_template_free", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature, batch_size)

    def predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict purification procedures for a given reaction SMILES string.
        """
        return self._predict("purification_protocols", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_purification_protocols", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature, batch_size)

    def predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict a product given a reactant SMILES string using a template-free approach.
        """
        return self._predict("forward_reaction", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict products given reactant SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_forward_reaction", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature, batch_size)

    def predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict retrosynthetic procedures for a given reactants and products reaction SMILES string.
        """
        return self._predict("procedures_given_reactants_products", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_procedures_given_reactants_products", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature, batch_size)

    def predict_top_k_retro_templated(self, input_data, model_version='latest', top_k=10):
        """
        Child function to predict retrosynthetic procedures for a given reactants and products reaction SMILES string.
        """
        return self._predict("top_k_retro_templated", input_data, 'smiles', model_version, top_k)

    def extract_reaction_procedure_jsons_from_text(self, input_data, model_version='latest', output_data_format='zip', upload_to_external_storage=True):
        """
        Extracts reaction procedure JSONs from a list of text passages.
    
        This function processes a list of text passages to extract reaction procedures,
        which are then returned as JSON objects. The function validates the input data
        and sends a POST request to a specified endpoint for processing.
    
        Parameters:
            input_data : list of str
                A list of text passages from which to extract reaction procedures.
                Each passage should be a string.
        
            model_version : str, optional
                The version of the model to use for extraction. Defaults to 'latest'.

            compress_input: bool, optional
                Whether to compress the input_data before sending post request. Defaults to True. 

            output_data_format: str, optional
                Format to return processed data. Must be one of XXXX

            upload_to_external_storage: bool, optional
                Whether to upload the output data to external storage. Defaults to True. 
    
        Returns:
            dict
                A dictionary containing the following keys:
                - 'input': The input data sent for processing.
                - 'output': The output data received from the processing endpoint.
                - 'status': The status of the processing request.
                - 'execution_time': The time taken to process the request in seconds.
    
        Raises:
            TypeError
                If 'input_data' is not a list of strings or if 'model_version' is not a string.
            ValueError
                If 'input_data' contains more than 300 passages or if any passage is longer than 8192 characters.
                If 'input_data_type' is not 'smiles' or 'rxn_smiles'.
        """
    
        if not isinstance(input_data, list) or not all(isinstance(item, str) for item in input_data):
            raise TypeError("The 'input_data' argument must be a list of strings.")
    
        if len(input_data) > 300:
            raise ValueError("The 'input_data' argument must not contain more than 300 passages.")
        
        for datum in input_data:
            if len(datum) > 8192:
                raise ValueError("Passages in the 'input_data' argument must not be longer than 8192 characters.") 
    
        if not isinstance(model_version, str):
            raise TypeError("The 'model_version' argument must be a string.")

        memory_file = io.BytesIO()
        with zipfile.ZipFile(memory_file, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
            zip_file.writestr('input_data.txt', input_data)
        memory_file.seek(0)
        zip_content = memory_file.getvalue()
        input_data = base64.b64encode(zip_content).decode('utf-8')
    
        post_request_data = {
            "endpoint": "reaction_extraction",
            "data": {
                "model_version": "extract_reaction_procedure_jsons_from_text",
                "input_data": input_data,
                "kwargs": {
                    "compress_input": True,
                    "output_data_format": 'zip',
                    "upload_to_external_storage": upload_to_external_storage
                }
            }
        }
    
        post_request_data = {"input": post_request_data}

        output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
    
        returned_data = {
            'input': post_request_data['input'],
            'job_id': output_data.json()['id'],
            'status': output_data.json()['stats']
        }
        
        return returned_data

    def _process_completed_response(response):
        """
        Process a completed API response, extract and decode the ZIP content,
        and return the processed JSON data.
    
        Args:
            response (requests.Response): The API response object.
    
        Returns:
            dict: The processed JSON data.
    
        Raises:
            ValueError: If the response status is not 'COMPLETED'.
        """
        if response.json()['status'] != 'COMPLETED':
            raise ValueError("Response status is not COMPLETED")
    
        encoded_zip = response.json()['output']
        zip_content = base64.b64decode(encoded_zip)
    
        with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as temp_zip:
            temp_zip.write(zip_content)
            temp_zip_path = temp_zip.name
    
        json_data = self._process_zip_file(temp_zip_path)
    
        Path(temp_zip_path).unlink()
        
        return list(json_data.values())
    
    def _process_zip_file(zip_path):
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
                    text = self._process_text_file(archive, file_name)
                    json_data = literal_eval(text)
                    combined_json_data[file_name] = json_data
        return combined_json_data
    
    def _process_text_file(archive, file_name):
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

    def check_status_extract_reaction_procedure_jsons_from_text(self, input_data, wait_to_complete = True, timeout = 1800):
            """
            Checks the status of reaction procedure extraction job.
        
            This function takes the the output data of an asynchronus request to extract reaction
            data and returns data based on user supplied inputs. 
        
            Parameters:
                input_data : XXX
                    The data returned from the asynchronus request to extract reaction data. 
    
                wait_to_complete : bool, optional
                    Whether to wait for the job to finish, and then return the extracted reaction data. 
                    Defaults to True. 
    
                timeout: float, optional
                    The amount of time to wait for the job to finish if wait_to_complete is True. Must be
                    a positive number less than 3600. Defaults to 1800. 
        
            Returns:
            
        
            Raises:
                
            """

            post_request_data = {
                "id": input_data["id"],
            }
    
            output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
        
            response_status = output_data.json()['status']
    
            if wait_to_complete == True:
                start = time.time()
                try:
                    while response_status not in ['COMPLETED', 'FAILED'] and time.time()-start < timeout:
                        time.sleep(5)
                        response = requests.post(status_url, 
                                        headers=headers)
                        response_status = response.json()['status']
                except KeyboardInterrupt:
                    print("\nInterrupted by user.")
    
                if response_status == 'COMPLETED':
                    output_data = self._process_completed_response(response.json())
                    
                    returned_data = {
                        'status': response_status,
                        'output': output_data
                    }

                else:
                    returned_data = {
                    'status': response_status,
                    'output': 'An error has occured'
                }
    
            else:
                returned_data = {
                    'status': response_status,
                }

                if response_status == 'COMPLETED':
                    output_data = self._process_completed_response(response.json())
                    returned_data['output'] = output_data
                    
            return returned_data
