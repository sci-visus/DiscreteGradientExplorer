# DiscreteGradientViewer

![Visualize construction of discrete gradient](sample.png)

## Overview

DiscreteGradientViewer is a tool designed to facilitate the visualization of Foreman's style of discrete Morse gradient fields. This tool enables users to visually inspect and analyze the construction of discrete gradient fields, aiding in understanding and exploring various discrete Morse theory concepts.

## Dependencies

Before using DiscreteGradientViewer, ensure you have the following dependencies installed:

- **freeglut**: A freely available alternative to the OpenGL Utility Toolkit (GLUT) library. It provides a simple API for creating windows, handling input, and more, making it suitable for graphics programming.
  
  Installation instructions:
  
  - **Linux** (Ubuntu/Debian):
    ```bash
    sudo apt-get install freeglut3-dev
    ```
    
  - **macOS** (via Homebrew):
    ```bash
    brew install freeglut
    ```
    
  - **Windows**: Binaries for Windows can be found online, or you can compile freeglut from source.

## Building

DiscreteGradientViewer uses CMake for building the project. Follow these steps to build DiscreteGradientViewer:

1. Clone the repository:
   ```bash
   git clone https://github.com/sci-visus/DiscreteGradientExplorer.git
   ```

2. Navigate to the project directory:
   ```bash
   cd DiscreteGradientViewer
   ```

3. Create a build directory:
   ```bash
   mkdir build
   cd build
   ```

4. Generate build files using CMake:
   ```bash
   cmake ..
   ```

5. Build the project:
   ```bash
   make
   ```

## Usage

After successfully building the project, you can run the `umsc` tool from the command line. Here's how to use it:

```bash
./umsc input_mesh.off [optional_scaling_value]
```

- `input_mesh.off`: Replace this with the path to your desired input mesh file in .off format. The z-coordinate of the mesh is treated as the scalar value.

- `[optional_scaling_value]`: An optional scaling value that you can provide to multiply the z-coordinate scalar values by. This helps in adjusting the visualization scale as needed.

Upon running the tool, you'll be presented with a window displaying the visualization of the discrete gradient fields. To understand how to navigate and interact with the visualization, please refer to the following files:

- **Instructions.txt**: This file provides a detailed guide on using mouse and keyboard commands to navigate the visualization and control various aspects of the display.

- **Sequence.txt**: Explore this file to understand how the discrete gradient construction unfolds within the context of a lower-stars filtration of the input mesh. It provides insights into the steps involved in the construction process and helps in visualizing the discrete gradient field's evolution.

Feel free to experiment with different commands and explore the visualization to gain a deeper understanding of the discrete gradient construction process. Adjust the optional scaling value as needed to optimize the visualization according to your preferences.

## Contributing

Contributions to DiscreteGradientViewer are welcome! If you encounter any issues or have ideas for improvements, feel free to open an issue or submit a pull request on GitHub.

## License

This project is licensed under the [BSD-3](LICENSE). Feel free to use, modify, and distribute this software as per the terms of the license.

