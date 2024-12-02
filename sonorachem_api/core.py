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
                                     json=data,
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
                ("Warning: Atom mapping found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
            
            if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                ("Warning: Isotopes found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetIsotope(0)
            
            Chem.SanitizeMol(mol)
            
            return True
        except Exception as e:
            (f"Error: {e}")
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
                    ("Warning: Atom mapping found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
    
            for mol in reaction.GetProducts():
                if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                    ("Warning: Atom mapping found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
            
            for mol in reaction.GetReactants():
                if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                    ("Warning: Isotopes found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetIsotope(0)
    
            for mol in reaction.GetProducts():
                if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                    ("Warning: Isotopes found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetIsotope(0)
            
            Chem.rdChemReactions.SanitizeRxn(reaction)
            
            if reaction.Validate()[1] != 0:
                return False
            
            return True
        except Exception as e:
            (f"Error: {e}")
            return False

    def _verify_kwargs(self, input_data, input_data_type='smiles', model_version='latest', kwargs={}, batched=False):
        """
        Verifies and validates the input arguments and keyword arguments for the model.
        
        Args:
            input_data (str or list): The input data to be processed. If `batched` is True, 
                this should be a list of strings. Otherwise, it should be a single string.
            input_data_type (str, optional): The type of input data. Must be one of 'smiles' 
                or 'rxn_smiles'. Defaults to 'smiles'.
            model_version (str, optional): The version of the model to be used. 
                Defaults to 'latest'.
            kwargs (dict, optional): Additional keyword arguments for the model. 
                Defaults to {}.
            batched (bool, optional): If True, the input data is expected to be a list 
                of strings. Defaults to False.
        
        Returns:
            None
        
        Raises:
            TypeError: 
                - If any of the input arguments or keyword arguments have incorrect types.
            ValueError: 
                - If any of the input arguments or keyword arguments have invalid values.
        """
    
        sampling_method = kwargs.get('sampling_method', 'greedy')
        seq_length = kwargs.get('seq_length', 256)
        beam_size = kwargs.get('beam_size', 5)
        temperature = kwargs.get('temperature', 0.3)
        top_k = kwargs.get('top_k', 16)
        batch_size = kwargs.get('batch_size', 16)

        if input_data_type not in ['smiles', 'rxn_smiles']:
            raise ValueError("Invalid 'input_data_type'. Must be 'smiles', 'rxn_smiles'.")
    
        if batched:
            if not isinstance(input_data, list) or not all(isinstance(item, str) for item in input_data):
                raise TypeError("The 'input_data' argument must be a list of strings.")
            
            if not isinstance(batch_size, int):
                raise TypeError("The 'batch_size' argument must be an integer.")
            if batch_size <= 0 or batch_size > 256:
                raise ValueError("The 'batch_size' argument must be greater than 0 and less than or equal to 256.")
        
            if input_data_type == 'smiles':
                for smiles in input_data:
                    valid_smiles = self.is_valid_smiles(smiles)
                    if not valid_smiles:
                        raise ValueError(f"The SMILES string '{smiles}' is not valid.")
                        
            elif input_data_type == 'rxn_smiles':
                for rxn_smiles in input_data:
                    valid_rxn_smiles = self.is_valid_reaction_smiles(rxn_smiles)
                    if not valid_rxn_smiles:
                        raise ValueError(f"The reaction SMILES string '{rxn_smiles}' is not valid.")
    
        else:  
            if not isinstance(input_data, str):
                raise TypeError("The 'input_data' argument must be a string.")
    
            if input_data_type == 'smiles':
                valid_smiles = self.is_valid_smiles(input_data)
                if not valid_smiles:
                    raise ValueError("The 'input_data' argument is not a valid SMILES string.")
                        
            elif input_data_type == 'rxn_smiles':
                valid_rxn_smiles = self.is_valid_reaction_smiles(input_data)
                if not valid_rxn_smiles:
                    raise ValueError("The 'input_data' argument is not a valid reaction SMILES string.")
    
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

        if not isinstance(temperature, (int, float)):
            raise TypeError("The 'temperature' argument must be an integer or a float.")
        if temperature <= 0:
            raise ValueError("The 'temperature' argument must be positive.")
    
        if not isinstance(top_k, int):
            raise TypeError("The 'top_k' argument must be an integer.")
        if top_k <= 0 or top_k > 128:
            raise ValueError("The 'top_k' argument must be greater than 0 and less than or equal to 128.")
                
    def _predict(self, endpoint, input_data, input_data_type='smiles', model_version='latest', kwargs={}):
        """
        Parent function to make a prediction.

        Args:
            endpoint (str): The endpoint to use for the prediction.
            input_data (str): input data.
            input_data_type (str, optional): The input data type used for input validation. 
                Must be one of 'smiles' or 'rxn_smiles'. Defaults to 'smiles'.
            model_version (str, optional): The version of the model to use. Defaults to 'latest'.
            kwargs (dict, optional): Arguments for predictions.

        Returns:
            dict: A dictionary containing the input, output, status, and execution time.

        Raises:
            ValueError: 
              - If arguments cannot be validated by _verify_kwargs.
        """

        self._verify_kwargs(input_data, input_data_type, model_version, kwargs, batched=False)
        
        post_request_data = {
            "endpoint": endpoint,
            "data": {
                "model_version": model_version,
                "batched": False,
                "input_data": input_data,
                "kwargs": kwargs
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

    def _batch_predict(self, endpoint, input_data,  input_data_type='smiles', model_version='latest', kwargs={}, batch_size=64):
        """
        Parent function to make batch predictions.

        Args:
            endpoint (str): The endpoint to use for the prediction.
            input_data (list of str): A list of SMILES strings.
            input_data_type (str, optional): The input data type used for input validation. 
                Must be one of 'smiles' or 'rxn_smiles'. Defaults to 'smiles'.
            model_version (str, optional): The version of the model to use. Defaults to 'latest'.
            kwargs (dict, optional): Arguments for predictions.
            batch_size (int, optional): The batch size for batch inference. Defaults to 64.

        Returns:
            dict: A dictionary containing the input, output, status, and execution time.

        Raises:
            ValueError: 
              - If arguments cannot be validated by _verify_kwargs.
        """

        kwargs["batch_size"] = batch_size

        self._verify_kwargs(input_data, input_data_type, model_version, kwargs, batched=True)
        
        post_request_data = {
            "endpoint": endpoint,
            "data": {
                "model_version": model_version,
                "input_data": input_data,
                "batched": True,
                "kwargs": kwargs
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
        return self._predict("procedures_retro_template_free", input_data, input_data_type='smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def batch_predict_procedures_retro_template_free(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict retrosynthetic procedures for SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_procedures_retro_template_free", input_data, input_data_type='smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict purification procedures for a given reaction SMILES string.
        """
        return self._predict("purification_protocols", input_data, input_data_type='rxn_smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def batch_predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_purification_protocols", input_data, input_data_type='rxn_smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict a product given a reactant SMILES string using a template-free approach.
        """
        return self._predict("forward_reaction", input_data, input_data_type='smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def batch_predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict products given reactant SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_forward_reaction", input_data, input_data_type='smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict retrosynthetic procedures for a given reactants and products reaction SMILES string.
        """
        return self._predict("procedures_given_reactants_products", input_data, input_data_type='rxn_smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def batch_predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3, batch_size=64):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_procedures_given_reactants_products", input_data, input_data_type='rxn_smiles', model_version=model_version, kwargs={'sampling_method': sampling_method, 'seq_length': seq_length, 'beam_size': beam_size, 'temperature': temperature})

    def predict_top_k_retro_templated(self, input_data, model_version='latest', top_k=16, rerank_by_reactants=True, use_saguarochem=True, use_custom_data=False):
        """
        Child function to predict top k most likely reactants from templates given a product.
        """
        return self._predict("top_k_retro_templated", input_data, input_data_type='smiles', model_version=model_version, kwargs={'top_k': top_k, 'rerank': rerank_by_reactants, 'use_saguarochem': use_saguarochem, 'use_custom_data': use_custom_data})

    def retrieve_top_k_similar_reactions_rxn_smiles(self, input_data, model_version='latest', top_k=16, use_saguarochem=True, use_custom_data=False):
        """
        Child function to retrieve top k most similar reactions given a reaction SMILES.
        """
        return self._predict("top_k_similar_reactions_rxn_smiles", input_data, input_data_type='rxn_smiles', model_version=model_version, kwargs={'top_k': top_k, 'use_saguarochem': use_saguarochem, 'use_custom_data': use_custom_data})

    def retrieve_top_k_similar_reactions_reactants(self, input_data, model_version='latest', top_k=16, use_saguarochem=True, use_custom_data=False):
        """
        Child function to retrieve top k most similar reactions given reactant SMILES.
        """
        return self._predict("top_k_similar_reactions_reactants", input_data, input_data_type='smiles', model_version=model_version, kwargs={'top_k': top_k, 'use_saguarochem': use_saguarochem, 'use_custom_data': use_custom_data})

    def retrieve_top_k_similar_reactions_products(self, input_data, model_version='latest', top_k=16, use_saguarochem=True, use_custom_data=False):
        """
        Child function to retrieve top k most similar reactions given a product SMILES.
        """
        return self._predict("top_k_similar_reactions_products", input_data, input_data_type='smiles', model_version=model_version, kwargs={'top_k': top_k, 'use_saguarochem': use_saguarochem, 'use_custom_data': use_custom_data})

    def extract_reaction_procedure_jsons_from_text(self, input_data, model_version='latest', compress_input=True, output_data_format='binary', upload_to_external_storage=True):
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
                Format to return processed data. Must be one of 'raw_output' or 'binary'. Defaults
                to 'binary'.

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

        if compress_input:            
            memory_file = io.BytesIO()
            with zipfile.ZipFile(memory_file, 'w', compression=zipfile.ZIP_DEFLATED) as zip_file:
                zip_file.writestr('input_data.txt', str(input_data))
            memory_file.seek(0)
            zip_content = memory_file.getvalue()
            input_data = base64.b64encode(zip_content).decode('utf-8')
    
        post_request_data = {
            "endpoint": "reaction_extraction",
            "data": {
                "model_version": model_version,
                "input_data": input_data,
                "kwargs": {
                    "compress_input": compress_input,
                    "output_data_format": output_data_format,
                    "upload_to_external_storage": upload_to_external_storage
                }
            }
        }
    
        post_request_data = {"input": post_request_data}

        output_data = self._send_post_request(self._base_url, self._headers, post_request_data)
    
        returned_data = {
            'input': post_request_data['input'],
            'job_id': output_data['id'],
            'status': output_data['status']
        }
        
        return returned_data

    def _process_completed_response(self, response, output_data_format='raw_output'):
        """
        Process a completed API response, extract and decode the ZIP content,
        and return the processed JSON data.
    
        Args:
            response (requests.Response): The API response object.
            output_data_format (str): The requested output data format.
    
        Returns:
            dict: The processed JSON data.
    
        Raises:
            ValueError: If the response status is not 'COMPLETED'.
        """
        if response['status'] != 'COMPLETED':
            raise ValueError("Response status is not COMPLETED")

        if output_data_format == 'raw_output':
            output_data = response['output']
            
        if output_data_format == 'binary':
            encoded_zip = response['output']
            zip_content = base64.b64decode(encoded_zip)
        
            with tempfile.NamedTemporaryFile(suffix='.zip', delete=False) as temp_zip:
                temp_zip.write(zip_content)
                temp_zip_path = temp_zip.name
        
            json_data = self._process_zip_file(temp_zip_path)
        
            Path(temp_zip_path).unlink()
            output_data = list(json_data.values())

        else:
            output_data = []
        
        return output_data
    
    def _process_zip_file(self, zip_path):
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
    
    def _process_text_file(self, archive, file_name):
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
                input_data : dict
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

            top_level_input = input_data["input"]['data']
            kwargs = top_level_input.get("kwargs", {})
            compress_input = kwargs.get("compress_input", False)
            output_data_format = kwargs.get("output_data_format", "binary")
            upload_to_external_storage = kwargs.get("upload_to_external_storage", False)
            job_id = input_data.get('job_id', None)

            post_request_data = {
                "endpoint": "check_reaction_extraction_status",
                "data": {
                    "model_version": "latest",
                    "input_data": job_id,
                    "kwargs": {}
                }
            }

            post_request_data = {"input": post_request_data}
    
            response = self._send_post_request(self._base_url, self._headers, post_request_data)
            response_status = response['status']
    
            if wait_to_complete == True:
                start = time.time()
                try:
                    while response_status not in ['COMPLETED', 'FAILED'] and time.time()-start < timeout:
                        time.sleep(5)
                        response = self._send_post_request(self._base_url, self._headers, post_request_data)
                        response_status = response['status']
                        
                except KeyboardInterrupt:
                    raise Exception("\nInterrupted by user.")
    
                if response_status == 'COMPLETED':
                    output_data = self._process_completed_response(response, output_data_format)
                    
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
                    output_data = self._process_completed_response(response, output_data_format)
                    returned_data['output'] = output_data
                    
            return returned_data
