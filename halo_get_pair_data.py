from yt.utilities.sdf import load_sdf
import math
import urllib2

HUBBLE_CONST = 0.688062

sim_file = '/home/b/Projects/DarkSky/Catalog/data_release/halos/ds14_a_halos_filter_1e15_1.0000'

#urllib2.urlopen('http://darksky.slac.stanford.edu/simulations/ds14_a/halos/ds14_a_halos_1.0000')

save_directory = '/home/b/Projects/DarkSky/Catalog/DSS_scripts/'

print "File Read.\n\n\n\n\n\n\n\n\nLeave Alone Please\n\n\n\n\n\n"

sdf_data = load_sdf(sim_file)

f_pairs_data = open(save_directory+'reduced_halo_pairs_full_data.txt','w')

#header describes the format of how a pair is stored
f_pairs_data.write('# pair_id\n')
f_pairs_data.write('# ax ay az avx avy avz amvir ar200b\n')
f_pairs_data.write('# bx by bz bvx bvy bvz bmvir br200b\n')

f_pairs_data.close()

f_pairs = open(save_directory+'reduced_halo_pairs.txt','r')
#f_pairs.next() #skip header line
i=0
for line in f_pairs:

    if i % 1000 == 0:
        print i

    halo_a = int(line.split()[0])
    halo_b = int(line.split()[1])

    ax = sdf_data['x'][halo_a]/HUBBLE_CONST
    ay = sdf_data['y'][halo_a]/HUBBLE_CONST
    az = sdf_data['z'][halo_a]/HUBBLE_CONST

    avx = sdf_data['vx'][halo_a]
    avy = sdf_data['vy'][halo_a]
    avz = sdf_data['vz'][halo_a]

    amvir = sdf_data['mvir'][halo_a]/HUBBLE_CONST
    ar200b = sdf_data['r200b'][halo_a]/HUBBLE_CONST

    bx = sdf_data['x'][halo_b]/HUBBLE_CONST
    by = sdf_data['y'][halo_b]/HUBBLE_CONST
    bz = sdf_data['z'][halo_b]/HUBBLE_CONST

    bvx = sdf_data['vx'][halo_b]
    bvy = sdf_data['vy'][halo_b]
    bvz = sdf_data['vz'][halo_b]

    bmvir = sdf_data['mvir'][halo_b]/HUBBLE_CONST
    br200b = sdf_data['r200b'][halo_b]/HUBBLE_CONST

    f_pairs_data = open(save_directory+'reduced_halo_pairs_full_data.txt','a')
    f_pairs_data.write(str(i)+'\n')
    f_pairs_data.write(str(halo_a)+' '+str(ax)+' '+str(ay)+' '+str(az)+' '+str(avx)+' '+str(avy)+' '+str(avz)+' '+str(amvir)+' '+str(ar200b)+' '+'\n')
    f_pairs_data.write(str(halo_b)+' '+str(bx)+' '+str(by)+' '+str(bz)+' '+str(bvx)+' '+str(bvy)+' '+str(bvz)+' '+str(bmvir)+' '+str(br200b)+' '+'\n')

    i+=1

    f_pairs_data.close()


f_pairs.close()
