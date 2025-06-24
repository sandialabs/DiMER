# script created using Sandia AIchat

import requests
import os

def download_file_from_zenodo(doi_or_id, output_directory):
    # Construct the Zenodo API URL
    url = f'https://zenodo.org/api/records/{doi_or_id}'
    
    # Make a GET request to the Zenodo API
    response = requests.get(url)
    
    if response.status_code == 200:
        # Parse the JSON response
        data = response.json()
        
        # Extract the files from the response
        files = data['files']
        
        # Create the output directory if it doesn't exist
        os.makedirs(output_directory, exist_ok=True)
        
        for file in files:
            file_url = file['links']['self']  # Get the file download link
            file_name = file['key']  # Get the file name
            
            # Download the file
            file_response = requests.get(file_url, stream=True)
            if file_response.status_code == 200:
                # Save the file to the specified directory
                with open(os.path.join(output_directory, file_name), 'wb') as f:
                    for chunk in file_response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f'Downloaded: {file_name}')
            else:
                print(f'Failed to download {file_name}: {file_response.status_code}')
    else:
        print(f'Failed to retrieve record: {response.status_code}')

# Example usage
doi_or_id = '15724972'  # Replace with your DOI or record ID
output_directory = './'  # Directory to save downloaded files
download_file_from_zenodo(doi_or_id, output_directory)
