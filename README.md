# DiscreteGradientViewer

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
   git clone https://github.com/your_username/DiscreteGradientViewer.git
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

After successfully building the project, you can run DiscreteGradientViewer from the command line. Here's how to use it:

```bash
./umsc <data.off> 1
```

Upon running the tool, you'll be presented with a window displaying the visualization of the mesh. Use mouse and keyboard inputs to navigate and interact with the visualization as needed. Refer to Instructions.txt, and Sequence.txt for keyboard and mouse commands. 

## Contributing

Contributions to DiscreteGradientViewer are welcome! If you encounter any issues or have ideas for improvements, feel free to open an issue or submit a pull request on GitHub.

## License

This project is licensed under the [BSD-3](LICENSE). Feel free to use, modify, and distribute this software as per the terms of the license.

