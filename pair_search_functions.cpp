#include "pair_search_functions.h"

double magnitude(cart_t cart){
  double mag;

  mag = sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);

  return mag;
}

//Projection of a onto b
cart_t projection(cart_t a, cart_t b){

  cart_t proj;

  proj.x = ( (a.x*b.x) + (a.y*b.y) + (a.z*b.z) ) / ( magnitude(b)*magnitude(b) ) * b.x;
  proj.y = ( (a.x*b.x) + (a.y*b.y) + (a.z*b.z) ) / ( magnitude(b)*magnitude(b) ) * b.y;
  proj.z = ( (a.x*b.x) + (a.y*b.y) + (a.z*b.z) ) / ( magnitude(b)*magnitude(b) ) * b.z;

  return proj;
}

cart_t sep_projection(cart_t a, cart_t b){

  cart_t p_sep, norm_comp;

  norm_comp = projection(a,b);

  p_sep.x = a.x - norm_comp.x;
  p_sep.y = a.y - norm_comp.y;
  p_sep.z = a.z - norm_comp.z;

  return p_sep;
}

cart_t sph_to_cart(sph_t sph){

  cart_t cart;

  cart.x = sin(sph.theta)*cos(sph.phi);
  cart.y = sin(sph.theta)*sin(sph.phi);
  cart.z = cos(sph.theta);

  return cart;
}

sph_t cart_to_sph(cart_t cart){

  sph_t sph;

  sph.theta = atan(sqrt(cart.x*cart.x + cart.y*cart.y)/cart.z);
  sph.phi = atan(cart.y/cart.x);

  return sph;
}

//relative velocity of a from b
cart_t get_rel_v(halo_t halo_a, halo_t halo_b){

  cart_t v;

  v.x = halo_a.vel.x - halo_b.vel.x;
  v.y = halo_a.vel.y - halo_b.vel.y;
  v.z = halo_a.vel.z - halo_b.vel.z;

  return v;
}

//relative position of a from b
cart_t get_rel_p(halo_t halo_a, halo_t halo_b){

  cart_t p;

  p.x = halo_a.pos.x - halo_b.pos.x;
  p.y = halo_a.pos.y - halo_b.pos.y;
  p.z = halo_a.pos.z - halo_b.pos.z;

  return p;
}

halo_t halo_t_parser(std::string str_input){
  int i;
  halo_t halo;
  std::string str_working[9];

  std::stringstream str_stream(str_input);

  if (str_stream.good()){
    for( i=0; i<N_HALO_ATTR; i++){
      str_stream >> str_working[i];
    }
  }

  halo.index = std::stoll(str_working[0].c_str());
  halo.pos.x = atof(str_working[1].c_str());
  halo.pos.y = atof(str_working[2].c_str());
  halo.pos.z = atof(str_working[3].c_str());
  halo.vel.x = atof(str_working[4].c_str());
  halo.vel.y = atof(str_working[5].c_str());
  halo.vel.z = atof(str_working[6].c_str());
  halo.mvir = atof(str_working[7].c_str());
  halo.r200b = atof(str_working[8].c_str());

  return halo;
}

bounds_t get_range_input(std::string type){

  bounds_t bound;
  double input, range;

  std::string units;

  if(type=="separation"){
    units = "(Mpc)";
  }

  if(type=="velocity"){
    units = "(km/s)";
  }

  if(type=="mass_a" || type=="mass_b"){
    units = "(Msun)";
  }

  std::cout << "Input projected " << type << " " << units <<": ";
  std::cin >> input;

  std::cout << "Input projected range " << type << " " << units <<": ";
  std::cin >> range;

  bound.up = input + range;
  bound.low = input - range;

  return bound;
}

double probability(std::string type, double mean, double sigma, double value){
	double prob;
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
void print_halo(halo_t halo){

  std::cout << halo.index << " " << halo.pos.x << " " << halo.pos.y << " " << halo.pos.z << " " << halo.vel.x << " " << halo.vel.y << " " << halo.vel.z << " " << halo.mvir << " " << halo.r200b  << std::endl;

  return;
}

void save_halo(halo_t halo, std::ofstream& data){

  data << halo.index << " " << halo.pos.x << " " << halo.pos.y << " " << halo.pos.z << " " << halo.vel.x << " " << halo.vel.y << " " << halo.vel.z << " " << halo.mvir << " " << halo.r200b  << std::endl;

  return;
}
