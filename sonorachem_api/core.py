import requests
import json
import time
from typing import Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem

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
            # Validate URL
            if not url.startswith(('http://', 'https://')):
                raise ValueError("Invalid URL. Must start with 'http://' or 'https://'")
    
            # Send the POST request
            response = requests.post(url, 
                                     headers=headers, 
                                     json=data,  # Use json parameter for automatic JSON encoding
                                     timeout=timeout, 
                                     verify=verify_ssl)
            
            # Raise an HTTPError for bad responses (4xx and 5xx status codes)
            response.raise_for_status()
    
            # Attempt to parse and return the JSON response
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
        # Prepare the input data with the API key
        input_data = {"api_key": self._api_key}
    
        # Send a POST request to the usage URL with the provided headers and input data
        output_data = self._send_post_request(self._usage_url, self._headers, input_data)
    
        # Return the output data containing the API usage statistics
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
            # Parse the SMILES string
            mol = Chem.MolFromSmiles(smiles)
            
            # Check if the molecule is valid
            if mol is None:
                return False
            
            # Check for atom mapping and remove if present
            if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                print("Warning: Atom mapping found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
            
            # Check for isotopes and remove if present
            if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                print("Warning: Isotopes found and removed.")
                for atom in mol.GetAtoms():
                    atom.SetIsotope(0)
            
            # Sanitize the molecule
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
            # Parse the reaction SMILES
            reaction = AllChem.ReactionFromSmarts(reaction_smiles)
            
            # Check if the reaction is valid
            if reaction is None:
                return False
            
            # Check for atom mapping and remove if present
            for mol in reaction.GetReactants() + reaction.GetProducts():
                if any(atom.GetAtomMapNum() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Atom mapping found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetAtomMapNum(0)
            
            # Check for isotopes and remove if present
            for mol in reaction.GetReactants() + reaction.GetProducts():
                if any(atom.GetIsotope() != 0 for atom in mol.GetAtoms()):
                    print("Warning: Isotopes found and removed.")
                    for atom in mol.GetAtoms():
                        atom.SetIsotope(0)
            
            # Sanitize the reaction
            Chem.SanitizeRxn(reaction)
            
            # Check if the reaction can be performed
            if reaction.Validate()[1] != 0:
                return False
            
            return True
        except Exception as e:
            print(f"Error: {e}")
            return False

    def predict_procedures_retro_template_free(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts retrosynthetic procedures for a given SMILES string using a template-free approach.
  
      This function takes a SMILES string and predicts potential retrosynthetic
      procedures. The function allows for different sampling methods and parameters to 
      control the prediction process.
  
      Args:
          input_data (str): A SMILES string.
          model_version (str, optional): The version of the model to use. Defaults to 'latest'.
          sampling_method (str, optional): The method used for sampling predictions. 
              Must be one of 'top_k', 'greedy', or 'sampling'. Defaults to 'top_k'.
          seq_length (int, optional): The maximum sequence length for the model input. 
              Defaults to 512.
          beam_size (int, optional): The beam size for beam search (if applicable). 
              Defaults to 5.
          temperature (float, optional): The temperature parameter for controlling randomness 
              in sampling. Must be a positive float. Higher values increase randomness. 
              Defaults to 1.0.
  
      Returns:
          list of dict: or something, fill in later when format is finalized
  
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
      # Validate input parameters
      if not isinstance(input_data, str):
          raise ValueError("The 'input_data' argument must be a string.")

      if not isinstance(model_version, str):
          raise ValueError("The 'model_version' argument must be a string.")
      
      if sampling_method not in ['top_k', 'greedy', 'sampling']:
          raise ValueError("Invalid sampling method. Must be 'top_k', 'greedy', or 'sampling'.")

      if not isinstance(seq_length, int):
          raise ValueError("seq_length must be an integer.")
      if seq_length <= 0 or seq_length > 512:
          raise ValueError("seq_length must be greater than 0 and less than or equal to 512.")
  
      if not isinstance(beam_size, int):
          raise ValueError("beam_size must be an integer.")
      if beam_size <= 0 or beam_size > 16:
          raise ValueError("beam_size must be greater than 0 and less than or equal to 16.")
        
      if beam_size == 1:
          sampling_method = 'greedy'
      
      if temperature <= 0:
          raise ValueError("Temperature must be a positive float.")
          
      valid_smiles = self.is_valid_smiles(input_data)
      if not valid_smiles:
            raise ValueError("The 'input_data' argument is not a valid SMILES string.")

      post_request_data = {
                    "endpoint": "procedures_retro_template_free",
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

      start = time.time()
      output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
      returned_data = {'input': post_request_data, 'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_purification_protocols(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts purification procedures for given reaction SMILES string.
  
      This function takes a reaction SMILES string and predicts potential purification 
      procedures. The function allows for different sampling methods and parameters to 
      control the prediction process.
  
      Args:
          input_data (str): A reaction SMILES string.
          model_version (str, optional): The version of the model to use. Defaults to 'latest'.
          sampling_method (str, optional): The method used for sampling predictions. 
              Must be one of 'top_k', 'greedy', or 'sampling'. Defaults to 'top_k'.
          seq_length (int, optional): The maximum sequence length for the model input. 
              Defaults to 512.
          beam_size (int, optional): The beam size for beam search (if applicable). 
              Defaults to 5.
          temperature (float, optional): The temperature parameter for controlling randomness 
              in sampling. Must be a positive float. Higher values increase randomness. 
              Defaults to 1.0.
  
      Returns:
          list of dict: or something, fill in later when format is finalized
  
      Raises:
          ValueError: 
            - If `input_data` is not a string.
            - If `model_version` is not a string.
            - If `sampling_method` is not one of 'top_k', 'greedy', or 'sampling'.
            - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
            - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
            - If `temperature` is not a positive float.
      """
      # Validate input parameters
      if not isinstance(input_data, str):
          raise ValueError("The 'input_data' argument must be a string.")

      if not isinstance(model_version, str):
          raise ValueError("The 'model_version' argument must be a string.")
      
      if sampling_method not in ['top_k', 'greedy', 'sampling']:
          raise ValueError("Invalid sampling method. Must be 'top_k', 'greedy', or 'sampling'.")

      if not isinstance(seq_length, int):
          raise ValueError("seq_length must be an integer.")
      if seq_length <= 0 or seq_length > 512:
          raise ValueError("seq_length must be greater than 0 and less than or equal to 512.")
  
      if not isinstance(beam_size, int):
          raise ValueError("beam_size must be an integer.")
      if beam_size <= 0 or beam_size > 16:
          raise ValueError("beam_size must be greater than 0 and less than or equal to 16.")
        
      if beam_size == 1:
          sampling_method = 'greedy'
      
      if temperature <= 0:
          raise ValueError("Temperature must be a positive float.")

      valid_reaction = self.is_valid_reaction_smiles(input_data)
      if not valid_reaction:
            raise ValueError("The 'input_data' argument is not a valid reaction SMILES string.")

      post_request_data = {
                    "endpoint": "purification_protocols",
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

      start = time.time()
      output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
      returned_data = {'input': post_request_data, 'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_forward_reaction(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts a product given reactant SMILES using a template-free approach.
  
      This function takes a reactant SMILES strings and predicts potential products
      The function allows for different sampling methods and parameters to control the 
      prediction process.
  
      Args:
          input_data (str): Reactant SMILES strings.
          model_version (str, optional): The version of the model to use. Defaults to 'latest'.
          sampling_method (str, optional): The method used for sampling predictions. 
              Must be one of 'top_k', 'greedy', or 'sampling'. Defaults to 'top_k'.
          seq_length (int, optional): The maximum sequence length for the model input. 
              Defaults to 512.
          beam_size (int, optional): The beam size for beam search (if applicable). 
              Defaults to 5.
          temperature (float, optional): The temperature parameter for controlling randomness 
              in sampling. Must be a positive float. Higher values increase randomness. 
              Defaults to 1.0.
  
      Returns:
          list of dict: or something, fill in later when format is finalized
  
      Raises:
          ValueError: 
            - If `input_data` is not a string.
            - If `model_version` is not a string.
            - If `sampling_method` is not one of 'top_k', 'greedy', or 'sampling'.
            - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
            - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
            - If `temperature` is not a positive float.
      """
      # Validate input parameters
      if not isinstance(input_data, str):
          raise ValueError("The 'input_data' argument must be a string.")

      if not isinstance(model_version, str):
          raise ValueError("The 'model_version' argument must be a string.")
      
      if sampling_method not in ['top_k', 'greedy', 'sampling']:
          raise ValueError("Invalid sampling method. Must be 'top_k', 'greedy', or 'sampling'.")

      if not isinstance(seq_length, int):
          raise ValueError("seq_length must be an integer.")
      if seq_length <= 0 or seq_length > 512:
          raise ValueError("seq_length must be greater than 0 and less than or equal to 512.")
  
      if not isinstance(beam_size, int):
          raise ValueError("beam_size must be an integer.")
      if beam_size <= 0 or beam_size > 16:
          raise ValueError("beam_size must be greater than 0 and less than or equal to 16.")
        
      if beam_size == 1:
          sampling_method = 'greedy'
      
      if temperature <= 0:
          raise ValueError("Temperature must be a positive float.")

      valid_smiles = self.is_valid_smiles(input_data)
      if not valid_smiles:
            raise ValueError("The 'input_data' argument is not a valid SMILES string.")

      post_request_data = {
                    "endpoint": "forward_reaction",
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

      start = time.time()
      output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
      returned_data = {'input': post_request_data, 'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_procedures_given_reactants_products(self, input_data, model_version='latest', sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts retrosynthetic procedures for given reactants and products SMILES string.
  
      This function takes a reaction SMILES string with reactants and products and predicts potential procedures.
      The function allows for different sampling methods and parameters to control the prediction process.
  
      Args:
          input_data (str): A reaction SMILES string.
          model_version (str, optional): The version of the model to use. Defaults to 'latest'.
          sampling_method (str, optional): The method used for sampling predictions. 
              Must be one of 'top_k', 'greedy', or 'sampling'. Defaults to 'top_k'.
          seq_length (int, optional): The maximum sequence length for the model input. 
              Defaults to 512.
          beam_size (int, optional): The beam size for beam search (if applicable). 
              Defaults to 5.
          temperature (float, optional): The temperature parameter for controlling randomness 
              in sampling. Must be a positive float. Higher values increase randomness. 
              Defaults to 1.0.
  
      Returns:
          list of dict: or something, fill in later when format is finalized
  
      Raises:
          ValueError: 
            - If `input_data` is not a string.
            - If `model_version` is not a string.
            - If `sampling_method` is not one of 'top_k', 'greedy', or 'sampling'.
            - If `seq_length` is not an integer, or if it is not greater than 0 or exceeds 512.
            - If `beam_size` is not an integer, or if it is not greater than 0 or exceeds 16.
            - If `temperature` is not a positive float.
      """
      # Validate input parameters
      if not isinstance(input_data, str):
          raise ValueError("The 'input_data' argument must be a string.")

      if not isinstance(model_version, str):
          raise ValueError("The 'model_version' argument must be a string.")
      
      if sampling_method not in ['top_k', 'greedy', 'sampling']:
          raise ValueError("Invalid sampling method. Must be 'top_k', 'greedy', or 'sampling'.")

      if not isinstance(seq_length, int):
          raise ValueError("seq_length must be an integer.")
      if seq_length <= 0 or seq_length > 512:
          raise ValueError("seq_length must be greater than 0 and less than or equal to 512.")
  
      if not isinstance(beam_size, int):
          raise ValueError("beam_size must be an integer.")
      if beam_size <= 0 or beam_size > 16:
          raise ValueError("beam_size must be greater than 0 and less than or equal to 16.")
        
      if beam_size == 1:
          sampling_method = 'greedy'
      
      if temperature <= 0:
          raise ValueError("Temperature must be a positive float.")

      invalid_reaction = self.is_valid_reaction_smiles(input_data)
      if not valid_reaction:
            raise ValueError("The 'input_data' argument is not a valid reaction SMILES string.")

      post_request_data = {
                    "endpoint": "procedures_given_reactants_products",
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

      start = time.time()
      output_data = self._send_post_request(self._runpod_url, self._headers, post_request_data)
      returned_data = {'input': post_request_data, 'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}
        
      return returned_data
