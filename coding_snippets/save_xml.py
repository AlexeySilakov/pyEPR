import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom
import base64
from collections import OrderedDict

def save_mixed_dict_xml(data_dict, filename):
    """
    Save a dictionary with mixed data types to an XML file.
    - NumPy arrays are encoded as base64 strings with metadata
    - Other data types (str, int, float, bool) are saved as text
    - Nested dictionaries are handled recursively
    """

    def create_element(parent, key, value):
        """Create XML element for a value"""
        element = ET.SubElement(parent, 'item')
        element.set('key', str(key))

        if isinstance(value, np.ndarray):
            # Handle numpy arrays
            array_element = ET.SubElement(element, 'numpy_array')
            array_element.set('dtype', str(value.dtype))
            array_element.set('shape', str(value.shape))

            # Encode array data as base64
            data_b64 = base64.b64encode(value.tobytes()).decode('utf-8')
            data_element = ET.SubElement(array_element, 'data')
            data_element.text = data_b64

        elif isinstance(value, dict):
            # Handle nested dictionaries
            dict_element = ET.SubElement(element, 'dict')
            for k, v in value.items():
                create_element(dict_element, k, v)

        elif isinstance(value, (list, tuple)):
            # Handle lists and tuples
            list_element = ET.SubElement(element, 'list')
            list_element.set('type', type(value).__name__)
            for item in value:
                create_element(list_element, 'item', item)

        else:
            # Handle primitive types
            primitive_element = ET.SubElement(element, 'primitive')
            primitive_element.set('type', type(value).__name__)
            primitive_element.text = str(value)

    # Create root element
    root = ET.Element('dictionary')

    # Process each item in the dictionary
    for key, value in data_dict.items():
        create_element(root, key, value)

    # Create pretty-printed XML
    rough_string = ET.tostring(root, encoding='unicode')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    # Write to file
    with open(filename, 'w') as f:
        f.write(pretty_xml)

def load_mixed_dict_xml(filename):
    """
    Load a dictionary from XML file, reconstructing original structure and types.
    """

    def parse_element(element):
        """Parse XML element and reconstruct value"""
        # Get the key from the element
        key = element.get('key')

        # Check what type of element this is
        numpy_array = element.find('numpy_array')
        dict_element = element.find('dict')
        list_element = element.find('list')
        primitive_element = element.find('primitive')

        if numpy_array is not None:
            # Reconstruct numpy array
            dtype = numpy_array.get('dtype')
            shape = eval(numpy_array.get('shape'))  # Safe for simple shapes
            data_b64 = numpy_array.find('data').text
            data = base64.b64decode(data_b64)
            array = np.frombuffer(data, dtype=dtype)
            return key, array.reshape(shape)

        elif dict_element is not None:
            # Reconstruct dictionary
            result = {}
            for item in dict_element.findall('item'):
                k, v = parse_element(item)
                if k is not None:
                    result[k] = v
            return key, result

        elif list_element is not None:
            # Reconstruct list/tuple
            items = []
            for item in list_element.findall('item'):
                _, v = parse_element(item)
                items.append(v)
            list_type = list_element.get('type')
            if list_type == 'tuple':
                return key, tuple(items)
            else:
                return key, items

        elif primitive_element is not None:
            # Reconstruct primitive type
            value_type = primitive_element.get('type')
            value_text = primitive_element.text

            if value_type == 'int':
                return key, int(value_text)
            elif value_type == 'float':
                return key, float(value_text)
            elif value_type == 'bool':
                return key, value_text.lower() == 'true'
            else:
                return key, value_text

        return key, None

    # Parse XML file
    tree = ET.parse(filename)
    root = tree.getroot()

    # Reconstruct dictionary
    result = {}
    for item in root.findall('item'):
        key, value = parse_element(item)
        if key is not None and value is not None:
            result[key] = value

    return result

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
        'tags': ['python', 'numpy', 'xml'],
        'coordinates': (10.5, 20.3)
    }

    print("Original data:")
    for k, v in test_data.items():
        print(f"  {k}: {type(v)} - {v}")

    # Save to XML
    save_mixed_dict_xml(test_data, 'test_data.xml')

    # Load from XML
    loaded_data = load_mixed_dict_xml('test_data.xml')

    print("\nLoaded data:")
    for k, v in loaded_data.items():
        print(f"  {k}: {type(v)} - {v}")

    # Verify they're equal
    # print("\nData integrity check:")
    # print(f"Names equal: {test_data['name'] == loaded_data['name']}")
    # print(f"Arrays equal: {np.array_equal(test_data['scores'], loaded_data['scores'])}")
    # print(f"Nested dict equal: {test_data['metadata'] == loaded_data['metadata']}")
    # print(f"Lists equal: {test_data['tags'] == loaded_data['tags']}")
    # print(f"Tuples equal: {test_data['coordinates'] == loaded_data['coordinates']}")

    # Display the XML file content
    print("\nXML file content:")
    with open('test_data.xml', 'r') as f:
        print(f.read())