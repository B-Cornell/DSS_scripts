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
