import numpy as np
import json
import base64
from collections import OrderedDict

def save_mixed_dict_json(data_dict, filename):
    """
    Save a dictionary with mixed data types to a JSON file.
    - NumPy arrays are encoded as base64 strings with metadata
    - Other data types (str, int, float, bool) are saved as-is
    - Nested dictionaries are handled recursively
    """

    def serialize_value(value):
        """Convert value to JSON-serializable format"""
        if isinstance(value, np.ndarray):
            # Convert numpy array to JSON-serializable format
            return {
                'type': 'numpy_array',
                'dtype': str(value.dtype),
                'shape': value.shape,
                'data': base64.b64encode(value.tobytes()).decode('utf-8')
            }
        elif isinstance(value, dict):
            # Recursively handle nested dictionaries
            return {
                'type': 'dict',
                'data': {k: serialize_value(v) for k, v in value.items()}
            }
        elif isinstance(value, (list, tuple)):
            # Handle lists and tuples
            return {
                'type': 'list',
                'data': [serialize_value(item) for item in value]
            }
        else:
            # For other types, return as-is
            return {
                'type': 'primitive',
                'data': value
            }

    # Serialize the entire dictionary
    serialized_data = {k: serialize_value(v) for k, v in data_dict.items()}

    # Save to JSON file
    with open(filename, 'w') as f:
        json.dump(serialized_data, f, indent=2)

def load_mixed_dict_json(filename):
    """
    Load a dictionary from JSON file, reconstructing original structure and types.
    """

    def deserialize_value(serialized_value):
        """Convert JSON data back to original types"""
        if serialized_value['type'] == 'numpy_array':
            # Reconstruct numpy array from base64 data
            dtype = np.dtype(serialized_value['dtype'])
            shape = serialized_value['shape']
            data = base64.b64decode(serialized_value['data'])
            array = np.frombuffer(data, dtype=dtype)
            return array.reshape(shape)
        elif serialized_value['type'] == 'dict':
            # Recursively handle nested dictionaries
            return {k: deserialize_value(v) for k, v in serialized_value['data'].items()}
        elif serialized_value['type'] == 'list':
            # Handle lists
            return [deserialize_value(item) for item in serialized_value['data']]
        else:
            # Return primitive value as-is
            return serialized_value['data']

    # Load from JSON file
    with open(filename, 'r') as f:
        serialized_data = json.load(f)

    # Deserialize back to original structure
    return {k: deserialize_value(v) for k, v in serialized_data.items()}

# Example usage and testing
if __name__ == "__main__":
    # Create test data
    test_data = {
        'name': 'John Doe',
        'age': 30,
        'height': 5.9,
        'is_student': False,
        'scores': np.array([85, 92, 78, 96], dtype=np.int32),
        'grades': np.array([[85, 92], [78, 96]], dtype=np.float64),
        'metadata': {
            'created': '2023-01-01',
            'version': 1.0,
            'active': True,
            'nested_array': np.array([1, 2, 3, 4, 5], dtype=np.int64)
        },
        'tags': ['python', 'numpy', 'json'],
        'coordinates': (10.5, 20.3)
    }

    print("Original data:")
    for k, v in test_data.items():
        print(f"  {k}: {type(v)} - {v}")

    # Save to JSON
    save_mixed_dict_json(test_data, 'test_data.json')

    # Load from JSON
    loaded_data = load_mixed_dict_json('test_data.json')

    print("\nLoaded data:")
    for k, v in loaded_data.items():
        print(f"  {k}: {type(v)} - {v}")

    # Verify they're equal
    print("\nData integrity check:")
    print(f"Names equal: {test_data['name'] == loaded_data['name']}")
    print(f"Arrays equal: {np.array_equal(test_data['scores'], loaded_data['scores'])}")
    print(f"Nested dict equal: {test_data['metadata'] } {loaded_data['metadata']}")
    print(f"Lists equal: {test_data['tags'] == loaded_data['tags']}")
    print(f"Tuples equal: {test_data['coordinates'] == loaded_data['coordinates']}")

    # Display the JSON file content
    print("\nJSON file content:")
    with open('test_data.json', 'r') as f:
        print(f.read())
