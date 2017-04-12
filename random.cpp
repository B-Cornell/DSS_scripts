/* There seems to be a problem with not canceling out orders of magnitude. If the mass is put in terms of E14 like it is supposed to, then the max probability is to the E-15 scale, but removing the E14 moves it up to E-1. I suspect that I don't have the variance correctly. I am going to have to ask about inputting the gaussian distribution. 

Else than that this works, I can just implement it in Jayke's code.
Now onto the alpha graph.
*/

#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cstdlib>
#include <vector>

double number;
double prob;
double PI = 3.1415926;
double largest;

double probability(std::string type, double mean, double sigma, double value){
	if (type == "tophat"){
		if (value<= mean + sigma and value>= mean - sigma){
			prob = 1;
		}
		else{
			prob = 0;
		}
	}
	if (type == "gaussian"){
	
		prob = (1/(sqrt(2*PI*sigma*sigma)))*exp(((mean-value)*(value-mean))/(2*sigma*sigma));
	}
 
	return prob;
}

int main(){
	
	std::cout << "Should be 1, 1, 0, lowish, high, low";
	std::cout <<std::endl;
	number = probability("tophat", 5, 10, 6);
	std::cout << number <<std::endl;
	number = probability("tophat", 1000, 2, 1000);
	std::cout << number <<std::endl;
	number = probability("tophat", 5, 10, 65);
	std::cout << number <<std::endl <<std::endl;


	number = probability("gaussian", 518.84599929, 108.195150286, 518.84599929);
	std::cout << "velocity    " << number <<std::endl;
	number = probability("gaussian", 3.55, 1.95, 3.55);
	largest = number;
	std::cout << "mass    " << number <<std::endl;
	number = probability("gaussian", 1, .3, 1);
	largest = number*largest;
	std::cout << "separation    " << number <<std::endl;
	largest = number*largest;
	std::cout <<std::endl << std::endl << largest << std::endl;
	return 0;
}
