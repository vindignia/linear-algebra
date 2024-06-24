# Description 
This package creates the matrix representation of spin quantum-mechanical operators and of some symmetry operators. Moreover, it finds the common eigenvectors of commuting hermitian matrices (operators).  


# Download this package (as artifact)

1. **Navigate to the Actions tab of this  GitHub repository**:
    - Click on the `Actions` tab at the top of the repository page.

2. **Select the (last) Workflow Run**:
    - Find the workflow run that contains the build you want to download. This is usually under the "All workflows" section.
    - Click on the specific workflow run.

3. **Download the package artifact**:
    - Scroll down to the `Artifacts` section in the workflow run summary.
    - Click on the "built-package" artifact to download it. This will download a ZIP file containing the built package (`.tar.gz` or `.whl` file).

## Install the package locally

Once you have downloaded the artifact, follow these steps to install the package:

1. **Extract the ZIP file**:
    - Unzip the downloaded file to access the `dist/` directory, which contains the `.tar.gz` or `.whl` files.

2. **Install the Package**:
    - Use `pip` to install the package. Run one of the following commands, depending on the file type:

For a `.whl` file:
```sh
pip install [PATH_TO_FILE]/linearalgebra-[VERSION]-py3-none-any.whl
```

For a `.tar.gz` file:
```sh
pip install [PATH_TO_FILE]/linearalgebra-[VERSION].tar.gz
```

By following these instructions, you can easily download and install and use this package in your local environment.
