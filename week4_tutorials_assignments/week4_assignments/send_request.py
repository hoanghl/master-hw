import requests 

def send_request(to_ig: bool=False, model_name: str="", ig_name: str="", user_type: str=None):
    """
    A helper function that sends some requests to a standalone inference service or an inference graph.
    Args:
        to_ig: If true, the request should be sent to an inference graph (ig), otherwise, directly to an inference service
        model_name: The name of the inference service if the request is sent to a separate inference service directly
        ig_name: The name of the inference graph if the request is sent to an inference graph first
        user_type: Include a "userType" field in the request if defined
    """
    request_data = {
        'parameters': {'content_type': 'pd'}, 
        'inputs': [{'name': 'season', 'shape': [2], 'datatype': 'UINT64', 'data': [1, 1]}, 
                   {'name': 'holiday', 'shape': [2], 'datatype': 'UINT64', 'data': [0, 0]}, 
                   {'name': 'workingday', 'shape': [2], 'datatype': 'UINT64', 'data': [0, 0]}, 
                   {'name': 'weather', 'shape': [2], 'datatype': 'UINT64', 'data': [1, 1]}, 
                   {'name': 'temp', 'shape': [2], 'datatype': 'FP64', 'data': [9.84, 9.02]}, 
                   {'name': 'atemp', 'shape': [2], 'datatype': 'FP64', 'data': [14.395, 13.635]}, 
                   {'name': 'humidity', 'shape': [2], 'datatype': 'UINT64', 'data': [81, 80]}, 
                   {'name': 'windspeed', 'shape': [2], 'datatype': 'FP64', 'data': [0.0, 0.0]}, 
                   {'name': 'hour', 'shape': [2], 'datatype': 'UINT64', 'data': [0, 1]}, 
                   {'name': 'day', 'shape': [2], 'datatype': 'UINT64', 'data': [1, 1]}, 
                   {'name': 'month', 'shape': [2], 'datatype': 'UINT64', 'data': [1, 1]}]
        }
    headers = {}
    kserve_gateway_url = "http://kserve-gateway.local:30200"
    if to_ig:
        headers["Host"] = f"{ig_name}.kserve-inference.example.com"
        if user_type != None:
            request_data["userType"] = user_type
        url = kserve_gateway_url
    else:
        headers["Host"] = f"{model_name}.kserve-inference.example.com"
        url = f"{kserve_gateway_url}/v2/models/{model_name}/infer"

    result = requests.post(url, json=request_data, headers=headers)
    print(result.json())

