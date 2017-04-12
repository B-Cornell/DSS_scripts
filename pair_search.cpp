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
  double azimuthal [ANGULAR_RES];
  double mass_a_a_prob; 
  double mass_a_b_prob; 
  double mass_b_a_prob; 
  double mass_b_b_prob; 
  double vel_prob; 
  double sep_prob; 
  double total_prob;


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
  f_halo_data.open("reduced_halo_pairs_full_data.txt");

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

  ofstream pair_out; //pair output
  pair_out.open(save_directory+"pair_output.txt");

  ofstream angle_out;
  angle_out.open(save_directory+"angle_out.txt");

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
    if(pair[k].id !=19022 && pair[k].id !=22172 && pair[k].id !=29604 && pair[k].id !=33600 && pair[k].id !=47215 && pair[k].id !=52167 && pair[k].id !=60438 && pair[k].id !=62786 && pair[k].id !=72103 && pair[k].id !=72362 && pair[k].id !=91533 && pair[k].id !=97408 && pair[k].id !=106445 && pair[k].id !=115912 && pair[k].id !=130910 && pair[k].id !=132036 && pair[k].id !=136944 && pair[k].id !=137521 && pair[k].id !=139327 && pair[k].id !=146295 && pair[k].id !=152876 && pair[k].id !=160395 && pair[k].id !=160611 && pair[k].id !=164999 && pair[k].id !=174842 && pair[k].id !=180739 && pair[k].id !=184029 && pair[k].id !=188011 && pair[k].id !=188819 && pair[k].id !=190248 && pair[k].id !=190722 && pair[k].id !=193113 && pair[k].id !=193380 && pair[k].id !=196802 && pair[k].id !=203020 && pair[k].id !=203033 && pair[k].id !=214808 && pair[k].id !=215392 && pair[k].id !=220797 && pair[k].id !=220997 && pair[k].id !=223161 && pair[k].id !=230861 && pair[k].id !=235991 && pair[k].id !=243196 && pair[k].id !=250492 && pair[k].id !=251696 && pair[k].id !=259711 && pair[k].id !=261923 && pair[k].id !=262991 && pair[k].id !=273782 && pair[k].id !=274464 && pair[k].id !=275773 && pair[k].id !=277873 && pair[k].id !=280505 && pair[k].id !=288477 && pair[k].id !=289086 && pair[k].id !=291178 && pair[k].id !=292849 && pair[k].id !=292918 && pair[k].id !=294428 && pair[k].id !=305389 && pair[k].id !=309429 && pair[k].id !=312686 && pair[k].id !=313153 && pair[k].id !=314794 && pair[k].id !=315909 && pair[k].id !=323305 && pair[k].id !=326843 && pair[k].id !=327152 && pair[k].id !=336909 && pair[k].id !=338822 && pair[k].id !=351828 && pair[k].id !=356319 && pair[k].id !=358086 && pair[k].id !=358636 && pair[k].id !=361307 && pair[k].id !=362745 && pair[k].id !=367027 && pair[k].id !=367104 && pair[k].id !=371068 && pair[k].id !=381905 && pair[k].id !=382799 && pair[k].id !=388192 && pair[k].id !=391007 && pair[k].id !=392033 && pair[k].id !=393385){
      calculation_pair = temp_halo(pair[k]);
    for(i = 0; sph.theta <= double(PI); i++){
      sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES)); // Range for theta is 0 to pi. 
      divisor = (ANGULAR_RES) * sin(sph.theta);
      total_divisor_points += divisor;
      points = int(divisor);
 
      sph.phi = 0.0;
    
      for( j = 0; sph.phi <= double(2*PI) and points>0; j++){
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
	  }
	  else {
	    azimuthal[99-i] += double(total_prob)*(divisor/(double(points)*double(points)));
	  }
	}
        }
      }
    
    calculation_pair.prob = area_counter / total_divisor_points;
    pair[k].prob = area_counter / total_divisor_points; //area that works divided by the surface area of the sphere
    area_counter = 0;
    
    

          
    //Print pair attributes
    cout << pair[k].id << endl;
    print_halo(calculation_pair.a);
    print_halo(calculation_pair.b);
    cout << "probability: " << pair[k].prob << endl;
    cout << "------------------------------------------" << endl;
       


} 
}         //outputting the angles to a file
          angle_out << "#" << endl;
          for( l = 0; l<ANGULAR_RES; l++){
            for( m = 0; m<ANGULAR_RES*2; m++){
              if(sphere[l][m] != '0'){
                angle_out << double(PI)/double(ANGULAR_RES) * l << " " << double(PI)/double(ANGULAR_RES) * m << endl;
              }
            }
          }

        
  }

  pair_out.close();
  angle_out.close();
  

  cout << "100%\nComplete."<< endl;

//##################################################################

  cout << endl << "Total Pairs: " << pair_count << endl << endl;
  cout << "(id)                  pair_id" << endl;
  cout << "(halo a attributes)   aindex ax ay az avx avy avz amvir ar200b" << endl;
  cout << "(halo b attributes)   bindex bx by bz bvx bvy bvz bmvir br200b" << endl;
  for(i = 0; i <51 ; i++){
    sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES));
    cout << azimuthal[i] << ",  ";
  }
  cout << endl ;

  return 0;
}
