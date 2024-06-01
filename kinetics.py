import numpy as np
import matplotlib.pyplot as plt

# 아레니우스 방정식을 이용하여 반응 속도 상수를 계산합니다
def arrhenius_equation(A, Ea, T):
    R = 8.314  # 기체 상수
    return A * np.exp(-Ea / (R * T))

class Chemical:
    def __init__(self, name, initial_concentration, coefficient, order):
        self.name = name
        self.concentration = initial_concentration
        self.coefficient = coefficient
        self.order = order
        self.delta = 0
        self.history = []

    def update_concentration(self, dt):
        self.concentration += self.delta * dt
        self.history.append(self.concentration)

    def apply_pressure_change(self, pressure_factor):
        self.concentration *= pressure_factor ** self.coefficient

class Reaction:
    def __init__(self, reactants, products, A_forward, Ea_forward, A_backward, Ea_backward, temperature):
        self.reactants = reactants
        self.products = products
        self.A_forward = A_forward
        self.Ea_forward = Ea_forward
        self.A_backward = A_backward
        self.Ea_backward = Ea_backward
        self.temperature = temperature
        self.isEqu = False
        self.equilibrium_timesteps = []
        self.equilibrium_concentrations = {c.name: [] for c in self.reactants + self.products}
        self.material_increase_steps = []
        self.temperature_change_steps = []
        self.pressure_change_steps = []

    def is_equilibrium(self, threshold=0.0001):
        for chemical in self.reactants + self.products:
            if abs(chemical.delta) > threshold:
                return False
        return True

    def update_reaction_rates(self):
        k_forward = arrhenius_equation(self.A_forward, self.Ea_forward, self.temperature)
        k_backward = arrhenius_equation(self.A_backward, self.Ea_backward, self.temperature)

        s_f = np.prod([chem.concentration ** chem.order for chem in self.reactants])
        s_b = np.prod([chem.concentration ** chem.order for chem in self.products])

        rate_forward = k_forward * s_f
        rate_backward = k_backward * s_b

        for chemical in self.reactants:
            chemical.delta = (-rate_forward + rate_backward) * chemical.coefficient
        for chemical in self.products:
            chemical.delta = (rate_forward - rate_backward) * chemical.coefficient

    def run_simulation(self, time_steps, dt, material_increase_steps, temperature_change_steps, pressure_change_steps):
        self.material_increase_steps = material_increase_steps
        self.temperature_change_steps = temperature_change_steps
        self.pressure_change_steps = pressure_change_steps
        
        for step in range(time_steps):
            for increase_step in material_increase_steps:
                if increase_step[0] == step:
                    for chem in self.reactants:
                        if chem.name == increase_step[1]:
                            chem.concentration += increase_step[2]

            for temp_change_step in temperature_change_steps:
                if temp_change_step[0] == step:
                    self.temperature += temp_change_step[1]

            for pressure_change_step in pressure_change_steps:
                if pressure_change_step[0] == step:
                    pressure_factor = pressure_change_step[1]
                    for chem in self.reactants + self.products:
                        chem.apply_pressure_change(pressure_factor)

            self.update_reaction_rates()

            if self.is_equilibrium() and not self.isEqu:
                self.equilibrium_timesteps.append(step)
                for chem in self.reactants + self.products:
                    self.equilibrium_concentrations[chem.name].append(chem.concentration)
                self.isEqu = True
            elif not self.is_equilibrium():
                self.isEqu = False

            for chem in self.reactants + self.products:
                chem.update_concentration(dt)
                
    def print_results(self):
        for k in range(len(self.material_increase_steps)):
            eqs = 0
            material_increase_step = self.material_increase_steps[k][0]
            for i in range(len(self.equilibrium_timesteps)):
                if material_increase_step < self.equilibrium_timesteps[i]:
                    eqs = self.equilibrium_timesteps[i]
                    break

            initial_equilibrium = {chem.name: chem.history[material_increase_step] for chem in self.reactants + self.products}
            final_equilibrium = {chem.name: chem.history[eqs] for chem in self.reactants + self.products}

            delta = {i: abs(initial_equilibrium[i] - final_equilibrium[i]) for i in initial_equilibrium}
            min_value = min(delta.values())

            print_result = f"Increase Material [{self.material_increase_steps[k][1]}] at time step [{self.material_increase_steps[k][0]}]\nDelta Materials : "
            for i in delta:
                print_result += f" {i}={delta[i]:.7f} "
            print_result += "\n(Ratio : "
            for i in delta:
                print_result += f" {i}={delta[i]/min_value:.7f} "
            print_result += ")\n"
            print(print_result)

        for k in range(len(self.temperature_change_steps)):
            eqs = 0
            temperature_increase_step = self.temperature_change_steps[k][0]
            for i in range(len(self.equilibrium_timesteps)):
                if temperature_increase_step < self.equilibrium_timesteps[i]:
                    eqs = self.equilibrium_timesteps[i]
                    break

            initial_equilibrium = {chem.name: chem.history[temperature_increase_step] for chem in self.reactants + self.products}
            final_equilibrium = {chem.name: chem.history[eqs] for chem in self.reactants + self.products}

            delta = {i: abs(initial_equilibrium[i] - final_equilibrium[i]) for i in initial_equilibrium}
            min_value = min(delta.values())

            print_result = f"Increase Temperature [{self.temperature_change_steps[k][1]}] at time step [{self.temperature_change_steps[k][0]}]\nDelta Materials : "
            for i in delta:
                print_result += f" {i}={delta[i]:.7f} "
            print_result += "\n(Ratio : "
            for i in delta:
                print_result += f" {i}={delta[i]/min_value:.7f} "
            print_result += ")\n"
            print(print_result)

        for k in range(len(self.pressure_change_steps)):
            eqs = 0
            pressure_change_step = self.pressure_change_steps[k][0]
            for i in range(len(self.equilibrium_timesteps)):
                if pressure_change_step < self.equilibrium_timesteps[i]:
                    eqs = self.equilibrium_timesteps[i]
                    break

            initial_equilibrium = {chem.name: chem.history[pressure_change_step] for chem in self.reactants + self.products}
            final_equilibrium = {chem.name: chem.history[eqs] for chem in self.reactants + self.products}

            delta = {i: abs(initial_equilibrium[i] - final_equilibrium[i]) for i in initial_equilibrium}
            min_value = min(delta.values())

            print_result = f"Change Pressure [Factor {self.pressure_change_steps[k][1]}] at time step [{self.pressure_change_steps[k][0]}]\nDelta Materials : "
            for i in delta:
                print_result += f" {i}={delta[i]:.7f} "
            print_result += "\n(Ratio : "
            for i in delta:
                print_result += f" {i}={delta[i]/min_value:.7f} "
            print_result += ")\n"
            print(print_result)

reactants = [
    Chemical('A', 1.0, 3, 1),
    Chemical('B', 2.3, 1, 1)
]

products = [
    Chemical('C', 0.0, 2, 1)
]

reaction = Reaction(
    reactants=reactants,
    products=products,
    A_forward=1e10,
    Ea_forward=50e3,
    A_backward=5e9,
    Ea_backward=45e3,
    temperature=273.15
)

time_steps = 1000
dt = 0.001
material_increase_steps = [[200, 'A', 0.5]]
temperature_change_steps = [[500, 35]]
pressure_change_steps = [[800, 0.8]]

reaction.run_simulation(time_steps, dt, material_increase_steps, temperature_change_steps, pressure_change_steps)
reaction.print_results()

# 결과 시각화
plt.figure(figsize=(10, 6))

for chem in reactants + products:
    plt.plot(chem.history, label=chem.name)

#for t in reaction.equilibrium_timesteps:
    #plt.axvline(x=t, color='g', linestyle='--', linewidth=0.7)
for t in reaction.material_increase_steps:
    plt.axvline(x=t[0], color='g', linestyle='--', linewidth=0.7)
for t in reaction.temperature_change_steps:
    plt.axvline(x=t[0], color='g', linestyle='--', linewidth=0.7)
for t in reaction.pressure_change_steps:
    plt.axvline(x=t[0], color='g', linestyle='--', linewidth=0.7)
plt.xlabel('Time Step')
plt.ylabel('Concentration')
#plt.rc('font', size=20)        # 기본 폰트 크기
#plt.rc('axes', labelsize=20)   # x,y축 label 폰트 크기
plt.rc('xtick', labelsize=10)  # x축 눈금 폰트 크기 
#plt.rc('ytick', labelsize=20)  # y축 눈금 폰트 크기
#plt.rc('legend', fontsize=20)  # 범례 폰트 크기
#plt.rc('figure', titlesize=50) # figure title 폰트 크기
plt.legend() 
plt.title('Le Chatelier\'s Principle Simulation')
plt.show()
