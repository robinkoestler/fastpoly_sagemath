{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "1. Install the NTL library using this [link](https://libntl.org/download.html). Use these [instructions](https://libntl.org/doc/tour-unix.html) to install NTL in /usr/local \n",
    "```bash\n",
    "   % gunzip ntl-11.5.1.tar.gz\n",
    "   % tar xf ntl-11.5.1.tar\n",
    "   % cd ntl-11.5.1/src\n",
    "   % ./configure \n",
    "   % make\n",
    "   % make check\n",
    "   % sudo make install\n",
    "```\n",
    "2. Check if a C++ compiler is available by typing in the terminal\n",
    "```bash\n",
    "   % clang --version\n",
    "```\n",
    "3. Move to the cloned git repository (of this) and execute\n",
    "```bash\n",
    "% g++ -g -O2 -std=c++11 -pthread -march=native testing_ntl.cpp -o testing_ntl -lntl -lgmp -lm && ./testing_ntl\n",
    "```\n",
    "It will create an executeable file called testing_ntl, which asks you for two numbers $a,b$ and will compute $(a+1)(b+1)$.\n",
    "   1. The tag -o stands for the output file and && concatenates two bash commands.\n",
    "   2. Moreover, ./testing_ntl will execute the compiled program."
   ]
  }
 ],
 "metadata": {
  "language_info": {
   "name": "python"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
