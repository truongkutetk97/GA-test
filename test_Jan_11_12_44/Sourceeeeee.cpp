#include <stdio.h>
#include <math.h>
extern "C"
{
	__declspec(dllexport)  void DisplayHelloFromDLL()
	{
		printf("Hello from DLL truong dep trai !\n");
		
	}
	__declspec(dllexport)  int Addd(int a, int b)
	 {
		return a + b;
	 }
	__declspec(dllexport)  float subb(int a, int b)
	 {
		 return a - b;
	 }
}