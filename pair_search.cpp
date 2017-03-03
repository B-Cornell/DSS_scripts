#include "pair_search_functions.h"

using namespace std;

int main(){

  vector<pair_t> pair(N_PAIRS); // This is our giant vector where we store all the halo pairs in the heap

  cart_t obs, obs_sep, obs_vel, rel_v, rel_p;
  sph_t sph;
  bounds_t b_sep, b_vel, b_mass_a, b_mass_b;
  double area_counter=0;
  double divisor;
  double total_divisor_points;
  int points;
  double azimuthal [ANGULAR_RES-1];
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
    if (k%5000 == 0){
      cout <<  double(k)/double(N_PAIRS)*100 << '%' << endl;
    }

    // Reset the sphere array
    for( i = 0; i<ANGULAR_RES; i++){
      for( j = 0; j<ANGULAR_RES*2; j++){
        sphere[i][j] = '0';
      }
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
    for(i = 0; sph.theta <= double(PI); i++){
      sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES)); // Range for theta is 0 to pi. 
      divisor = (ANGULAR_RES) * sin(sph.theta);
      total_divisor_points += divisor;
      points = int(divisor);
 
      sph.phi = 0.0;
    
      for( j = 0; sph.phi <= double(2*PI) and points>0; j++){
        sph.phi = sph.phi + (double(PI*2.0))/points; // Range for phi is 0 to 2pi
   
        obs = sph_to_cart(sph); // Convert spherical coordinates to cartesian
	
        rel_p = get_rel_p(pair[k].a,pair[k].b); // Calculate relative position
        rel_v = get_rel_v(pair[k].a,pair[k].b); // Calculate relative velocity

        obs_vel = projection(rel_v,obs); // Calculate observed velocity
        obs_sep = sep_projection(rel_p,obs); // Calculate observed separation

        mass_a_a_prob = probability("gaussian", b_mass_a.up - ((b_mass_a.up - b_mass_a.low)/2.), ((b_mass_a.up - b_mass_a.low)/2.), pair[k].a.mvir);//probability that the mass of the first halo in the pair matches the first halo in the musketball

	mass_a_b_prob = probability("gaussian", b_mass_a.up - ((b_mass_a.up - b_mass_a.low)/2.), ((b_mass_a.up - b_mass_a.low)/2.), pair[k].b.mvir);//probability that the mass of the first halo in the pair matches the second halo in the musketball
	
	mass_b_a_prob = probability("gaussian", b_mass_b.up - ((b_mass_b.up - b_mass_b.low)/2.), ((b_mass_b.up - b_mass_b.low)/2.), pair[k].a.mvir);//probability that the mass of the second halo in the pair matches the first halo in the musketball

	mass_b_b_prob = probability("gaussian", b_mass_b.up - ((b_mass_b.up - b_mass_b.low)/2.), ((b_mass_b.up - b_mass_b.low)/2.), pair[k].b.mvir);//probability that the mass of the second halo in the pair matches the second halo in the musketball

        vel_prob = probability("gaussian", b_vel.up - ((b_vel.up - b_vel.low)/2.), ((b_vel.up - b_vel.low)/2.), magnitude(obs_vel));

	sep_prob = probability("gaussian", b_sep.up - ((b_sep.up - b_sep.low)/2.), ((b_sep.up - b_sep.low)/2.), magnitude(obs_sep));
	
	total_prob = (mass_a_a_prob*mass_b_b_prob+mass_a_b_prob*mass_b_a_prob)*vel_prob*sep_prob;



        sphere[i][j] = ' '; // Mark where on the sphere the criterion is fulfilled
        area_counter += total_prob*(divisor/(double(points)*double(points)));
	azimuthal[i] += total_prob*(divisor/(double(points)*double(points))); 
        }
      }
    }

    pair[k].prob = area_counter / total_divisor_points; //area that works divided by the surface area of the sphere
    area_counter = 0;

    //checks if there is atleast one angle that works
    for( i = 0; i<ANGULAR_RES; i++){	

      for( j = 0; j<ANGULAR_RES*2; j++){
        // Check if there is at least one angle for which the pair can be an analog
	
        if(sphere[i][j] != '0'){

          pair_count++;
          //Print pair attributes
          cout << pair[k].id << endl;
          print_halo(pair[k].a);
          print_halo(pair[k].b);
          cout << "probability: " << pair[k].prob << endl;
          cout << "------------------------------------------" << endl;

          pair_out << pair[k].id << endl;
          save_halo(pair[k].a,pair_out);
          save_halo(pair[k].b,pair_out);
          pair_out << pair[k].prob << endl; //store data in output file

          i = ANGULAR_RES;
          j = ANGULAR_RES*2;

  /*        //Print out the array
          for( i = 0; i<ANGULAR_RES; i++){
            for( j = 0; j<ANGULAR_RES*2; j++){
              cout << sphere[i][j];
            }
            cout << endl;
          }

*/
          //outputting the angles to a file
          angle_out << "#" << endl;
          for( l = 0; l<ANGULAR_RES; l++){
            for( m = 0; m<ANGULAR_RES*2; m++){
              if(sphere[l][m] != '0'){
                angle_out << double(PI)/double(ANGULAR_RES) * l << " " << double(PI)/double(ANGULAR_RES) * m << endl;
              }
            }
          }

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
  for(i = 0; sph.theta < double(PI); i++){
    sph.theta = sph.theta + (double(PI)/double(ANGULAR_RES));
    cout << azimuthal[i] << ",  ";
  }
  cout << endl;

  return 0;
}
