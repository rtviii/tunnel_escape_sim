from pprint import pprint
import numpy as np
import pandas as pd
import numpy as np
from scipy import signal

# open csv with pandas:
df = pd.read_csv('tunnel_data.csv', header=None)
x_tunnel = df[0].to_numpy()
r_tunnel = df[1].to_numpy()


def interp1(x, y, x_new):
    """
    1-D interpolation function similar to MATLAB's interp1.
    
    Parameters:
        x (array-like): The x-coordinates of the data points.
        y (array-like): The y-coordinates of the data points.
        x_new (array-like): The x-coordinates at which to evaluate the interpolated values.
        
    Returns:
        y_new (array): The interpolated values at x_new.
    """
    # Find indices of the nearest neighbors in x for each value in x_new
    indices = np.searchsorted(x, x_new)
    
    # Clip indices to ensure they are within the valid range
    indices = np.clip(indices, 1, len(x)-1)
    
    # Calculate the fractions indicating the position of x_new between the nearest neighbors
    fractions = (x_new - x[indices-1]) / (x[indices] - x[indices-1])
    
    # Perform linear interpolation
    y_new = y[indices-1] + fractions * (y[indices] - y[indices-1])
    
    return y_new


# #!----------------------------
# x_tunnel = [1,2,3,4,5]
# r_tunnel = [1,2,3,4,5]
# #!----------------------------

L = max(x_tunnel)
new_x    = np.linspace(min(x_tunnel), max(x_tunnel), len(x_tunnel))
new_r = interp1(x_tunnel, r_tunnel, new_x)

x_tunnel = new_x;
y_tunnel = new_r;

print("Savitzky interpolation")
smoothed_y = signal.savgol_filter(y_tunnel, window_length=55, polyorder=3, mode="nearest")

R_radius_scaling = 1
x_init = x_tunnel[10];
K = 1000; # % Number of simulation runs
D = 0.5;
dt=0.01; # %delta_t
# M = len(R);
M = 1
T_Iterations = 1000;
MFPT_MF_passage_time = np.zeros(( M,1 ));

# % Initialize variables
# % h = animatedline;
FPT_f_passage_time=np.zeros(( K,M ));
for n in range(0,M):
    print("Simulating for value of M: {}".format(M))
    for k in range(1,K):
        if k%100==0:
                print(k)
        x = x_init;
        y = 0.1;  
        t = 0;   
        # // % Define tunnel boundaries
        y_lower = lambda x: -interp1(x_tunnel, R_radius_scaling*y_tunnel, x);
        y_upper = lambda x: interp1(x_tunnel, R_radius_scaling*y_tunnel, x);   #* ?

        while x > 0 and x < L and t<T_Iterations:  
            dx = np.sqrt(2*D*dt)*np.random.randn()
            dy = np.sqrt(2*D*dt)*np.random.randn()

            x = x + dx;
            y = y + dy;

            if x <= 0:
                x = x - dx;

            if y < y_lower(x) or y > y_upper(x):
                y = y - dy;
            t = t + 1;
        if x>=L:
            FPT_f_passage_time[k,n] = t;

for n in range(1,M):
    FPT_f_passage_time = FPT_f_passage_time[:,n];
    FPT_f_passage_time = FPT_f_passage_time[FPT_f_passage_time != 0];

pprint(FPT_f_passage_time)
# * ?
# MFPT[n] = np.mean(FPT)

# signal.savgol_filter()
# #!==============================================================================

# L = max(x_tunnel)

# new_x = np.linspace(min(x_tunnel), max(x_tunnel), len(x_tunnel));   # new x values
# new_r = np.interp1(x_tunnel,r, new_x);   # interpolate new y values  % interpolate new r values