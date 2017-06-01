#include "pair_search_functions.h"

using namespace std;

int main(){

  vector<pair_t> pair(N_PAIRS); // This is our giant vector where we store all the halo pairs in the heap

  cart_t obs, obs_sep, obs_vel, rel_v, rel_p;
  sph_t sph;
  bounds_t b_sep, b_vel, b_mass_a, b_mass_b;
  pair_t calculation_pair;
  double area_counter=0;
  double divisor;
  double total_divisor_points;
  int points;
  double azimuthal [ANGULAR_RES], azitemp[ANGULAR_RES];
  double mass_a_a_prob; 
  double mass_a_b_prob; 
  double mass_b_a_prob; 
  double mass_b_b_prob; 
  double vel_prob; 
  double sep_prob; 
  double angle_prob;
  double total_prob;
  int rel_vel_angle;
  double z_vel, y_vel;
  int returning[2];
  double rel_vel[91];
  

  string halo_a_str, halo_b_str, pair_id_str, temp;

  int i,j,k,l,m, pair_count=0;

  string save_directory = "/home/Projects/Darksky/Catalog/DSS_Scripts/";

/*how many analogs

 * sphere[theta][phi]
 * Theta has a range of 0 to pi.
 * Phi has a range of 0 to 2pi.
 * Since phi covers twice the interval, we do ANGULAR_RES*2.
 * This array stores all the viewing angles that fulfil the search criterion.
 * Printing this array is the same as taking the surface of a sphere and flattening and stretching it into a squre.
 * Angular resolution determines the number of "pixels" on this sphere.
 */
  char sphere[ANGULAR_RES][ANGULAR_RES*2];

  b_sep = get_range_input("separation"); // Units: Mpc
  b_vel = get_range_input("velocity"); // Units: km/s
  b_mass_a = get_range_input("mass_a"); // Units: Msun (Solar Masses)
  b_mass_b = get_range_input("mass_b"); // Units: Msun (Solar Masses)

  cout << "------------------------------------------" << endl;

  ifstream f_halo_data;
  f_halo_data.open("reduced_halo_pairs_full_data_3500Kpc.txt");

  if (f_halo_data.is_open()){

    // Skip header
    for(i = 0; i < N_HEADER_LINES; i++){
      getline(f_halo_data,temp);
    }

    // This loop sticks everything into a giant vector of type pair
    for(i=0; i< N_PAIRS; i++){

      getline(f_halo_data,pair_id_str); //halo pair id
      getline(f_halo_data,halo_a_str); // halo a data
      getline(f_halo_data,halo_b_str); // halo b data

      pair[i].id = atoi(pair_id_str.c_str());
      pair[i].a = halo_t_parser(halo_a_str); // Parses pair.a data into halo_t retainer
      pair[i].b = halo_t_parser(halo_b_str); // Parses pair.b data into halo_t retainer

      if (i%10000 == 0){
        cout << "Processing... " <<  double(i)/double(N_PAIRS)*100 << '%' << endl;
      }
    }
    cout << "Processing... 100%\nComplete."<< endl;
  }

  else {
    cout << "Error: Cannot open file." << endl;
    return 1;
  }

  cout << "Searching for matching pairs." << endl;

  
  std::ofstream pyfile;
  pyfile.open("pyfile.csv");
  pyfile << "#HaloID, Mass 1, Mass 2, Separation, Velocity on axis, Velocity off axis, Returning, Relative Velocity Angle, Angle probs, Total prob";
  

  // Iterates over all the pairs
  for(k=0; k < N_PAIRS; k++){

    // Print progress in percentage
    if (k%1000 == 0){
      cout <<  double(k)/double(N_PAIRS)*100 << '%' << endl;
    }

    // Integrating over the sphere
    // The difference between steps in theta are the same as steps in phi
    sph.theta = 0.;
    total_divisor_points = 0.;
    
    // Mass check
    if( ( ((pair[k].a.mvir > b_mass_a.low) && (pair[k].a.mvir <  b_mass_a.up))    &&
          ((pair[k].b.mvir > b_mass_b.low) && (pair[k].b.mvir <  b_mass_b.up)) )  ||
        ( ((pair[k].a.mvir > b_mass_b.low) && (pair[k].a.mvir <  b_mass_b.up))    &&
          ((pair[k].b.mvir > b_mass_a.low) && (pair[k].b.mvir <  b_mass_a.up)) )  ){
    pair_count++;
    calculation_pair = temp_halo(pair[k]);
      
      
    pyfile << pair[k].id << ", " << calculation_pair.a.mvir << ", " << calculation_pair.b.mvir << ", " << calculation_pair.b.pos.z << ", " << calculation_pair.b.vel.z << ", " << calculation_pair.b.vel.y << ", ";   
      /*find the velocity directions*/
    if (calculation_pair.b.vel.z < 0.0){
      z_vel = -calculation_pair.b.vel.z;
      returning[0]++;
      pyfile << "1, ";
    }
    else {
      z_vel = calculation_pair.b.vel.z;
      returning[1]++;
      pyfile << "0, ";
      }
    if (calculation_pair.b.vel.y < 0.0){
      y_vel = -calculation_pair.b.vel.y;
    }
    else {
      y_vel = calculation_pair.b.vel.y;

    }
    rel_vel_angle = atan(y_vel/z_vel) * (180 / PI);
    rel_vel[rel_vel_angle]++;
      
    pyfile << rel_vel_angle<< ", " ;
      
      
    for(i = 0; sph.theta <= double(PI); i++){
      sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES)); // Range for theta is 0 to pi. 
      divisor = (ANGULAR_RES) * sin(sph.theta);
      total_divisor_points += divisor;
      points = int(divisor);
 
      sph.phi = 0.0;
    
      for( j = 0; sph.phi < double(2*PI) and points>0; j++){
        sph.phi = sph.phi + (double(PI*2.0))/points; // Range for phi is 0 to 2pi
   
        obs = sph_to_cart(sph); // Convert spherical coordinates to cartesian
	
        rel_p = get_rel_p(calculation_pair.a,calculation_pair.b); // Calculate relative position
        rel_v = get_rel_v(calculation_pair.a,calculation_pair.b); // Calculate relative velocity

        obs_vel = projection(rel_v,obs); // Calculate observed velocity
        obs_sep = sep_projection(rel_p,obs); // Calculate observed separation

        mass_a_a_prob = probability("gaussian", b_mass_a.up - ((b_mass_a.up - b_mass_a.low)/2.), ((b_mass_a.up - b_mass_a.low)/2.), calculation_pair.a.mvir);//probability that the mass of the first halo in the pair matches the first halo in the musketball

	mass_a_b_prob = probability("gaussian", b_mass_a.up - ((b_mass_a.up - b_mass_a.low)/2.), ((b_mass_a.up - b_mass_a.low)/2.), calculation_pair.b.mvir);//probability that the mass of the first halo in the pair matches the second halo in the musketball
	
	mass_b_a_prob = probability("gaussian", b_mass_b.up - ((b_mass_b.up - b_mass_b.low)/2.), ((b_mass_b.up - b_mass_b.low)/2.), calculation_pair.a.mvir);//probability that the mass of the second halo in the pair matches the first halo in the musketball

	mass_b_b_prob = probability("gaussian", b_mass_b.up - ((b_mass_b.up - b_mass_b.low)/2.), ((b_mass_b.up - b_mass_b.low)/2.), calculation_pair.b.mvir);//probability that the mass of the second halo in the pair matches the second halo in the musketball

        vel_prob = probability("gaussian", b_vel.up - ((b_vel.up - b_vel.low)/2.), ((b_vel.up - b_vel.low)/2.), magnitude(obs_vel));

	sep_prob = probability("gaussian", b_sep.up - ((b_sep.up - b_sep.low)/2.), ((b_sep.up - b_sep.low)/2.), magnitude(obs_sep));
	
	total_prob = (mass_a_a_prob*mass_b_b_prob+mass_a_b_prob*mass_b_a_prob)*vel_prob*sep_prob;

        
        area_counter += double(total_prob)*(double(divisor)/(double(points)*double(points)));
	if (azimuthal[i] == azimuthal[i]){
	  if (i <50){
	    azimuthal[i] += double(total_prob)*(divisor/(double(points)*double(points)));
	    
	    azitemp[i] += double(total_prob)*(divisor/(double(points)*double(points)));
	  }
	  else {
	    azimuthal[99-i] += double(total_prob)*(divisor/(double(points)*double(points)));
	    azitemp[99-i] += double(total_prob)*(divisor/(double(points)*double(points)));
	  }
	}
	
        }
      }
    
    calculation_pair.prob = area_counter / total_divisor_points;
    pair[k].prob = area_counter / total_divisor_points; //area that works divided by the surface area of the sphere
    area_counter = 0;
    i =0;
    for(i = 0; i <50 ; i++){
    pyfile << azitemp[i] << ",  ";
    azitemp[i] = 0;
    }
    pyfile << calculation_pair.prob << endl;

          
    //Print pair attributes
    cout << pair[k].id << endl;
    print_halo(calculation_pair.a);
    print_halo(calculation_pair.b);
    cout << "probability: " << pair[k].prob << endl;
    cout << "relative velocity angle: " << rel_vel_angle << endl;
    cout << "------------------------------------------" << endl;
       


} 
}
         

  pyfile.close();

  cout << "100%\nComplete."<< endl;

//##################################################################

  cout << endl << "Total Pairs: " << pair_count << endl << endl;
  cout << "(id)                  pair_id" << endl;
  cout << "(halo a attributes)   aindex ax ay az avx avy avz amvir ar200b" << endl;
  cout << "(halo b attributes)   bindex bx by bz bvx bvy bvz bmvir br200b" << endl;
  for(i = 0; i <50 ; i++){
    sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES));
    cout << azimuthal[i] << ",  ";
  }
  cout << endl << endl;
  for(i = 0; i < 91; i++){
    cout << rel_vel[i] << ",  ";
  
  }
  cout << endl << endl;
  cout << returning[0] << ",  " << returning[1] << endl;
  return 0;
}
