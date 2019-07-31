from fractions import gcd
import os

get_dir = "/Volumes/parash-bkp/cer2/"
local_dir = "~/Desktop/analyze/"

for i in range(5):
    for j in range(5):
        print("ucer {}, ffa16 {}".format(i, j))
        multiplier = min(gcd(i, 4), gcd(j, 4))
        a = int(i / multiplier)
        b = int(4 / multiplier)
        c = int(j / multiplier)

        filename = "{}-{}-{}-{}-{}-2".format(a, b-a, b, c, b-c)

        curr_dir = os.path.join(get_dir, filename)
        traj_file = os.path.join(curr_dir, "simulation/traj.dcd")
        start_file = os.path.join(curr_dir, "simulation/start.hoomdxml")

        curr_loc_dir = os.path.join(local_dir, "{}-{}".format(i, j))
        os.system("mkdir -p {}".format(curr_loc_dir))
        os.system("scp percus:{} {}".format(start_file, curr_loc_dir))
        os.system("scp percus:{} {}".format(traj_file, curr_loc_dir))
