import requests
import json
import time
from typing import Dict, Any, Optional

class SaguaroChemAPIWrapper:
    """
    Python wrapper for SaguaroChem API.
    """

    def __init__(
        self,
        api_key: str,
        base_url: Optional[str] = 'https://hwsvflxtqh.execute-api.us-east-2.amazonaws.com/denovochem_api_stage/runpod',
    ):
        """
        SaguaroChemAPIWrapper constructor.

        Args:
            api_key (str): an API key to access the service.
            base_url (str, optional): base url for the service. If not provided it will default to
                the De Novo Chem AWS server

        Examples:
            Initialize the wrapper by simply providing an API key:

            >>> from xxx import SaguaroChemAPIWrapper
            >>> saguarochem_api_wrapper = SaguaroChemAPIWrapper(api_key=api_key)
        """
        self._api_key = api_key
        self._base_url = base_url
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

    def set_api_key(self, api_key: str):
        """
        Set the API key.

        This method also rebuilds the headers.

        Args:
            api_key (str): an API key to access the service.

        Examples:
            Set an API key:

            >>> saguarochem_api_wrapper.set_api_key('API_KEY')
        """
        logger.info("Set API key to {}".format(api_key))
        self._api_key = api_key
        self._headers = self._construct_headers()
        
    # def check_api_usage():

    def predict_procedures_retro_template_free(self, smiles, sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts retrosynthetic procedures for given SMILES strings using a template-free approach.
  
      This function takes a list of SMILES strings and predicts potential retrosynthetic
      procedures. The function allows for different sampling methods and parameters to 
      control the prediction process.
  
      Args:
          smiles (list of str): A list of SMILES strings representing chemical compounds.
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
          list of list of str: A list of predicted retrosynthetic procedures for each input SMILES.
              Each procedure is represented as a list of strings describing synthetic steps.
  
      Raises:
          ValueError: If an invalid sampling_method is provided or if temperature is not positive.
      """
      # Validate input parameters
      if not isinstance(smiles, list) or not all(isinstance(s, str) for s in smiles):
          raise ValueError("The 'smiles' argument must be a list of strings.")
      
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

      # input_data = {
      #               "endpoint": "procedures_retro_template_free",
      #               "data": {
      #                       "inference_strategy": "N/A",
      #                       "smiles": smiles,
      #                       "kwargs": {
      #                                 "sampling_method": sampling_method,
      #                                 "seq_length": seq_length,
      #                                 "beam_size": beam_size,
      #                                 "temperature": temperature
      #                                 }
      #                       }
      #               }

      input_data = {
                    "input": {
                        "smiles": smiles[0]
                    }
                }

      start = time.time()
      output_data = self._send_post_request(self._base_url, self._headers, input_data)
      returned_data = {'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_purification_protocols(self, smiles, sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts purification procedures for given reaction SMILES strings.
  
      This function takes a list of reaction SMILES strings and predicts potential 
      purification procedures. The function allows for different sampling methods and 
      parameters to control the prediction process.
  
      Args:
          smiles (list of str): A list of SMILES strings representing chemical compounds.
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
          list of list of str: A list of predicted retrosynthetic procedures for each input SMILES.
              Each procedure is represented as a list of strings describing synthetic steps.
  
      Raises:
          ValueError: If an invalid sampling_method is provided or if temperature is not positive.
      """
      # Validate input parameters
      if not isinstance(smiles, list) or not all(isinstance(s, str) for s in smiles):
          raise ValueError("The 'smiles' argument must be a list of strings.")
      
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

      input_data = {
                    "endpoint": "purification_protocols",
                    "data": {
                            "inference_strategy": "N/A",
                            "smiles": smiles,
                            "kwargs": {
                                      "sampling_method": sampling_method,
                                      "seq_length": seq_length,
                                      "beam_size": beam_size,
                                      "temperature": temperature
                                      }
                            }
                    }

      start = time.time()
      output_data = self._send_post_request(self._base_url, self._headers, input_data)
      returned_data = {'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_forward_reaction(self, smiles, sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts products given reactants SMILES using a template-free approach.
  
      This function takes a list of SMILES strings and predicts potential products
      The function allows for different sampling methods and parameters to control the 
      prediction process.
  
      Args:
          smiles (list of str): A list of SMILES strings representing chemical compounds.
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
          list of list of str: A list of predicted retrosynthetic procedures for each input SMILES.
              Each procedure is represented as a list of strings describing synthetic steps.
  
      Raises:
          ValueError: If an invalid sampling_method is provided or if temperature is not positive.
      """
      # Validate input parameters
      if not isinstance(smiles, list) or not all(isinstance(s, str) for s in smiles):
          raise ValueError("The 'smiles' argument must be a list of strings.")
      
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

      input_data = {
                    "endpoint": "forward_reaction",
                    "data": {
                            "inference_strategy": "N/A",
                            "smiles": smiles,
                            "kwargs": {
                                      "sampling_method": sampling_method,
                                      "seq_length": seq_length,
                                      "beam_size": beam_size,
                                      "temperature": temperature
                                      }
                            }
                    }

      start = time.time()
      output_data = self._send_post_request(self._base_url, self._headers, input_data)
      returned_data = {'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
      
    def predict_procedures_given_reactants_products(self, smiles, sampling_method='greedy', seq_length=256, beam_size=5, temperature=0.3):
      """
      Predicts retrosynthetic procedures for given reactants and products SMILES strings.
  
      This function takes a list of reactants and products SMILES strings and predicts potential
      procedures. The function allows for different sampling methods and parameters to control 
      the prediction process.
  
      Args:
          smiles (list of str): A list of SMILES strings representing chemical compounds.
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
          list of list of str: A list of predicted retrosynthetic procedures for each input SMILES.
              Each procedure is represented as a list of strings describing synthetic steps.
  
      Raises:
          ValueError: If an invalid sampling_method is provided or if temperature is not positive.
      """
      # Validate input parameters
      if not isinstance(smiles, list) or not all(isinstance(s, str) for s in smiles):
          raise ValueError("The 'smiles' argument must be a list of strings.")
      
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

      input_data = {
                    "endpoint": "procedures_given_reactants_products",
                    "data": {
                            "inference_strategy": "N/A",
                            "smiles": smiles,
                            "kwargs": {
                                      "sampling_method": sampling_method,
                                      "seq_length": seq_length,
                                      "beam_size": beam_size,
                                      "temperature": temperature
                                      }
                            }
                    }

      start = time.time()
      output_data = self._send_post_request(self._base_url, self._headers, input_data)
      returned_data = {'output': output_data['output'], 'status': output_data['status'], 'execution_time': time.time()-start}

      return returned_data
