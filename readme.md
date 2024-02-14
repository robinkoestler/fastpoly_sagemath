1. Install the NTL library using this [link](https://libntl.org/download.html). Use these [instructions](https://libntl.org/doc/tour-unix.html) to install NTL in /usr/local
   ```bash
      gunzip ntl-11.5.1.tar.gz
      tar xf ntl-11.5.1.tar
      cd ntl-11.5.1/src
      ./configure 
      make
      make check
      sudo make install
   ```
2. Check if a C++ compiler is available by typing in the terminal.
   ```bash
      clang --version
   ```
3. Move to the cloned git repository (of this) and execute
   ```bash
      g++ -g -O2 -std=c++11 -pthread -march=native testing_ntl.cpp -o testing_ntl -lntl -lgmp -lm && ./testing_ntl
   ```
It will create an executable file called testing_ntl, which asks you for two numbers $a,b$ and will compute $(a+1)(b+1)$.
   1. The tag -o stands for the output file and && concatenates two bash commands.
   2. Moreover, ./testing_ntl will execute the compiled program.
4. (UNRELATED) Include C++ code in SageMath. [Link](https://doc.sagemath.org/html/en/thematic_tutorials/cython_interface.html#calling-code-from-a-compiled-library)
