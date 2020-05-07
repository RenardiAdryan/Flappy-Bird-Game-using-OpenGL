#include <iostream>
#include <stdio.h>
#include <cstdlib>

using namespace std;

int main()
{
	system("g++ GlobalVariable.h serial.cpp FlappyBird.cpp -lGL -lGLU -lglut -o flappybird");
}