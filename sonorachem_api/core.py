import requests
import json
import time
from typing import Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*')                                                                                                                                                           

class SonoraChemAPIWrapper:
    """
    Python wrapper for SonoraChem API.
    """

    def __init__(
        self,
        api_key: str,
        base_url: Optional[str] = 'https://hwsvflxtqh.execute-api.us-east-2.amazonaws.com/denovochem_api_stage',
    ):
        """
        SonoraChemAPIWrapper constructor.

        Args:
            api_key (str): an API key to access the service.
            base_url (str, optional): base url for the service. If not provided it will default to
                the De Novo Chem AWS server

        Examples:
            Initialize the wrapper by simply providing an API key:

            >>> from sonorachem_api import SonoraChemAPIWrapper
            >>> sonorachem_api_wrapper = SonoraChemAPIWrapper(api_key=api_key)
        """
        self._api_key = api_key
        self._base_url = base_url
        self._runpod_url = self._base_url + '/runpod'
        self._usage_url = self._base_url + '/get_usage'
        self._headers = self._construct_headers()

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
        self._runpod_url = self._base_url + '/runpod'
        self._usage_url = self._base_url + '/get_usage'

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
        input_data = {"api_key": self._api_key}
    
        output_data = self._send_post_request(self._usage_url, self._headers, input_data)
    
        return output_data['body']

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
        output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
        returned_data = {
            'input': post_request_data,
            'output': output_data['output'],
            'status': output_data['status'],
            'execution_time': time.time() - start
        }
    
        return returned_data

    def _batch_predict(self, endpoint, input_data, input_data_type='smiles', model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
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

        Returns:
            dict: A dictionary containing the input, output, status, and execution time.

        Raises:
            ValueError: 
              - If `input_data` is not a list of strings.
              - If `model_version` is not a string.
              - If `sampling_method` is not one of 'greedy', or 'sampling'.
              - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
              - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
              - If `temperature` is not a positive float.
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
                    "temperature": temperature
                }
            }
        }

        post_request_data = {"input": post_request_data}

        start = time.time()
        output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
        returned_data = {
            'input': post_request_data,
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

    def batch_predict_procedures_retro_template_free(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to batch predict retrosynthetic procedures for SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_procedures_retro_template_free", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict purification procedures for a given reaction SMILES string.
        """
        return self._predict("purification_protocols", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_purification_protocols", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict a product given a reactant SMILES string using a template-free approach.
        """
        return self._predict("forward_reaction", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to batch predict products given reactant SMILES strings using a template-free approach.
        """
        return self._batch_predict("batch_forward_reaction", input_data, 'smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to predict retrosynthetic procedures for a given reactants and products reaction SMILES string.
        """
        return self._predict("procedures_given_reactants_products", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)

    def batch_predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
        """
        Child function to batch predict purification procedures for reaction SMILES strings.
        """
        return self._batch_predict("batch_procedures_given_reactants_products", input_data, 'rxn_smiles', model_version, sampling_method, seq_length, beam_size, temperature)
