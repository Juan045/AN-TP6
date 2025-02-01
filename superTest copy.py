import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys

class IceMeltingSimulator:
    def __init__(self, nx: int, ny: int, T_ice: float):
        # Physical parameters
        self.L = 334000.0  # Latent heat of fusion [J/kg]
        self.cp_ice = 2090.0  # Specific heat capacity ice [J/kg·K]
        self.cp_water = 4186.0  # Specific heat capacity water [J/kg·K]
        self.k_ice = 2.22  # Thermal conductivity ice [W/m·K]
        self.k_water = 0.556  # Thermal conductivity water [W/m·K]
        self.rho = 1000.0  # Density [kg/m³]
        self.T_melt = 0.0  # Melting temperature [°C]
        self.T_water = 20.0  # Water temperature [°C]

        # Grid parameters
        self.nx, self.ny = nx, ny
        self.dx = 0.01 / (nx-1)
        self.dy = 0.01 / (ny-1)

        # Initialize arrays with explicit float dtype
        self.T = np.full((ny, nx), T_ice, dtype=np.float64)
        self.phase = np.zeros((ny, nx), dtype=np.float64)
        self.domain_mask = self.create_t_shape()

        # Set water regions
        self.T[~self.domain_mask] = self.T_water
        self.phase[~self.domain_mask] = 1.0

        # Print initial temperature distribution
        plt.figure(figsize=(8, 8))
        plt.imshow(self.T, cmap='coolwarm')
        plt.colorbar(label='Temperature (°C)')
        plt.title(f'Initial Temperature Distribution\n{nx}x{nx} grid, T_ice = {T_ice}°C')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.show()

        # Also print the numerical values
        print("\nTemperature Matrix:")
        print(self.T)
        # Precompute constants
        alpha_max = max(self.k_water/(self.rho * self.cp_water), 
                       self.k_ice/(self.rho * self.cp_ice))
        self.dt = 0.2 * min(self.dx, self.dy)**2 / alpha_max

        # Precompute thermal properties arrays
        self.k = np.full((ny, nx), self.k_ice, dtype=np.float64)
        self.k[~self.domain_mask] = self.k_water
        self.cp = np.full((ny, nx), self.cp_ice, dtype=np.float64)
        self.cp[~self.domain_mask] = self.cp_water

    def create_t_shape(self):
        mask = np.zeros((self.ny, self.nx), dtype=bool)
        vert_width = self.nx // 5
        center = self.nx // 2
        mask[:3*self.ny//4, center-vert_width//2:center+vert_width//2+1] = True
        mask[3*self.ny//4:, :] = True
        return mask

    def solve_timestep(self):
        # Compute Laplacian terms
        dx2 = np.zeros_like(self.T)
        dy2 = np.zeros_like(self.T)

        dx2[1:-1, 1:-1] = (self.T[1:-1, 2:] - 2*self.T[1:-1, 1:-1] + self.T[1:-1, :-2])/self.dx**2
        dy2[1:-1, 1:-1] = (self.T[2:, 1:-1] - 2*self.T[1:-1, 1:-1] + self.T[:-2, 1:-1])/self.dy**2

        # Update temperature
        dT = self.dt * self.k/(self.rho * self.cp) * (dx2 + dy2)
        self.T[self.domain_mask] += dT[self.domain_mask]

        # Update phase for melting points
        melting_mask = (self.T > self.T_melt) & (self.phase < 1.0) & self.domain_mask
        self.phase[melting_mask] = np.minimum(1.0, 
            self.phase[melting_mask] + self.dt * self.k_ice * (self.T[melting_mask] - self.T_melt) / 
            (self.L * self.rho * self.dx**2))

        # Update thermal properties
        self.k[self.phase >= 1.0] = self.k_water
        self.cp[self.phase >= 1.0] = self.cp_water

        # Calculate and return liquid fraction
        return np.sum(self.phase[self.domain_mask]) / np.sum(self.domain_mask)

    def run_until_half_melted(self, max_time=1000.0):
        time = 0.0
        check_interval = 100  # Check progress every N steps
        pbar = tqdm(total=100, desc="Simulating")
        last_fraction = 0

        while time < max_time:
            for _ in range(check_interval):
                liquid_fraction = self.solve_timestep()
                time += self.dt

                if liquid_fraction >= 0.5:
                    pbar.close()
                    return time

            # Update progress bar
            progress = int(100 * liquid_fraction)
            if progress > last_fraction:
                pbar.update(progress - last_fraction)
                last_fraction = progress

        pbar.close()
        return max_time

def run_comparison():
    grid_sizes = [50, 100, 200]
    ice_temps = [-10, -20, -30]
    results = {}

    for nx in grid_sizes:
        for T_ice in ice_temps:
            print(f"\nRunning simulation: {nx}x{nx} grid, T_ice = {T_ice}°C")
            sim = IceMeltingSimulator(nx, ny=nx, T_ice=T_ice)
            melt_time = sim.run_until_half_melted()
            results[(nx, T_ice)] = melt_time
            print(f"Completed in {melt_time:.2f} seconds")

    return results

# Run and display results
results = run_comparison()
print("\nResults:")
print("=" * 50)
for (nx, T_ice), melt_time in sorted(results.items()):
    print(f"Grid size: {nx}x{nx}, Initial temp: {T_ice}°C, Melt time: {melt_time:.2f}s")

# Plot results
plt.figure(figsize=(10, 6))
for T_ice in [-10, -20, -30]:
    grid_sizes = [k[0] for k in results.keys() if k[1] == T_ice]
    melt_times = [results[k] for k in results.keys() if k[1] == T_ice]
    plt.plot(grid_sizes, melt_times, 'o-', label=f'T_ice = {T_ice}°C')

plt.xlabel('Grid Size')
plt.ylabel('Time to 50% Melting (s)')
plt.title('Melting Time vs Grid Size for Different Initial Temperatures')
plt.legend()
plt.grid(True)
plt.show()
