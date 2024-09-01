import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import re

# Set the print options for NumPy arrays to show 3 decimal places
np.set_printoptions(precision=3, suppress=True)

# Step 1: Get user input
beam_length = float(input("Enter the total length of the beam (in meters): "))
num_elements = int(input("Enter the number of elements you want: "))
num_nodes = num_elements + 1

# Determine node positions
element_length = beam_length / num_elements
nodes = np.linspace(0, beam_length, num_nodes)

# Get support types
supports = []
intermediate_support = input(f"Are there any intermediate support(s)?(Yes/No): ").strip().lower()
if (intermediate_support=="yes"):
    for i in range(num_nodes):
        support_type = input(f"Enter the type of support at node {i+1} (Fixed, Roller, Hinge, None): ").strip().lower()
        supports.append(support_type)
    # print(supports)
elif (intermediate_support=="no"):
    for i in range(num_nodes):
        if (i==0):
            support_type = input(f"Enter the type of support at node {i+1} (Fixed, Roller, Hinge, None): ").strip().lower()
            supports.append(support_type)
        elif(i<=num_nodes-2):
            supports.append("none")
        elif(i==(num_nodes-1)):
            support_type = input(f"Enter the type of support at node {i+1} (Fixed, Roller, Hinge, None): ").strip().lower()
            supports.append(support_type)

# print(supports)

# Get number of loads
num_loads = int(input("Enter the number of loads: "))
loads = []

for _ in range(num_loads):
    load_type = input("Enter the type of load (point/UDL/moment): ").strip().lower()
    
    if load_type == 'point':
        load_value = float(input("Enter the value of the point load (in kN): "))
        load_direction = input("Enter the direction of the load (up/down): ").strip().lower()
        load_position = float(input("Enter the position of the point load from the start of the beam: "))
        loads.append({'type': 'point', 'value': load_value, 'direction': load_direction, 'position': load_position})
    
    elif load_type == 'udl':
        load_value = float(input("Enter the value of the UDL (in kN/m): "))
        load_start = float(input("Enter the start position of the UDL: "))
        load_end = float(input("Enter the end position of the UDL: "))
        load_direction = input("Enter the direction of the UDL (up/down): ").strip().lower()
        loads.append({'type': 'udl', 'value': load_value, 'start': load_start, 'end': load_end, 'direction': load_direction})
    
    elif load_type == 'moment':
        moment_value = float(input("Enter the value of the moment (in kNm): "))
        moment_direction = input("Enter the direction of the moment (clockwise/anticlockwise): ").strip().lower()
        moment_position = float(input("Enter the position of the moment from the start of the beam: "))
        moment_value = -moment_value if moment_direction == 'clockwise' else moment_value
        loads.append({'type': 'moment', 'value': moment_value, 'position': moment_position,'direction':moment_direction})

# Display nodes
print(f"Nodes (equal division): {np.round(nodes, 3)}")

# Step 2: Calculate element lengths
element_lengths = np.diff(nodes)
print(f"Element Lengths: {np.round(element_lengths, 3)}")

# Step 3: Generate stiffness matrices for each element
EI = 1  # Placeholder for EI, results will be in terms of EI
stiffness_matrices = []

for i, length in enumerate(element_lengths):
    k = (EI / length**3) * np.array([
        [12, 6*length, -12, 6*length],
        [6*length, 4*length**2, -6*length, 2*length**2],
        [-12, -6*length, 12, -6*length],
        [6*length, 2*length**2, -6*length, 4*length**2]
    ])
    stiffness_matrices.append(k)
    # Print each element's stiffness matrix
    print(f"Stiffness Matrix for Element {i+1}:")
    print(np.array2string(np.round(k, 3), separator=', '))  # Improved matrix printing with 3 decimals and converted to string for separation by comma

# Step 4: Assemble the global stiffness matrix
global_matrix_size = 2 * num_nodes #Each node has 2 degrees of freedom assuming beam is inextensible

global_stiffness_matrix = np.zeros((global_matrix_size, global_matrix_size))

for i in range(num_elements):
    k = stiffness_matrices[i]
    indices = [2*i, 2*i+1, 2*(i+1), 2*(i+1)+1]
    for row in range(4):
        for col in range(4):
            global_stiffness_matrix[indices[row], indices[col]] += k[row, col]

print("Global Stiffness Matrix:")
print(np.array2string(np.round(global_stiffness_matrix, 3), separator=', '))  # Improved matrix printing with 3 decimals

# Step 5: Define the displacement vector [u] and force vector [F]
displacement_vector = np.zeros(global_matrix_size, dtype=object) #object type because it has to hold a mix of string and symbolic data
force_vector = np.zeros(global_matrix_size, dtype=object)

# Apply boundary conditions and forces
for i, support in enumerate(supports):
    if support == 'fixed':
        displacement_vector[2*i] = 0  # Vertical displacement is zero
        displacement_vector[2*i+1] = 0  # Rotation is zero
        force_vector[2*i] = f'R{i+1}'  # Reaction force
        force_vector[2*i+1] = f'M{i+1}'  # Moment
    elif support == 'roller' or support == 'hinge':
        displacement_vector[2*i] = 0  # Vertical displacement is zero
        displacement_vector[2*i+1] = f'Theta{i+1}'  # Unknown rotation
        force_vector[2*i] = f'R{i+1}'  # Reaction force
        force_vector[2*i+1] = 0  # Moment is zero
    elif support == 'none':
        displacement_vector[2*i] = f'V{i+1}'
        displacement_vector[2*i+1] = f'Theta{i+1}'
        force_vector[2*i] = 0  # No reaction force
        force_vector[2*i+1] = 0  # No moment

# Calculate and accumulate the force vectors for each element
element_force_vectors = []

for i in range(num_elements):
    element_force_vector = np.zeros(global_matrix_size)
    
    # Apply loads to the element
    for load in loads:
        if load['type'] == 'point':
            if nodes[i] <= load['position'] < nodes[i+1]:
                L_left = load['position'] - nodes[i]
                L_right = nodes[i+1] - load['position']
                length = nodes[i+1] - nodes[i]

                F_left = abs((load['value'] * L_right**2) / (length**3 * (3 * L_left + L_right)))
                F_right = abs((load['value'] * L_left**2) / (length**3 * (3 * L_right + L_left)))

                M_left = abs((load['value'] * L_left * L_right**2) / (length**2))
                M_right = abs((load['value'] * L_right * L_left**2) / (length**2))

                if load['direction'] == 'down':
                    F_left = -F_left
                    F_right = -F_right
                    M_left=-M_left
                elif load['direction']=='up':
                    M_right=-M_right

                element_force_vector[2*i] += F_left
                element_force_vector[2*i+1] += M_left
                element_force_vector[2*(i+1)] += F_right
                element_force_vector[2*(i+1)+1] += M_right

        elif load['type'] == 'udl':
            # UDL load characteristics
            w = load['value']
            load_start = load['start']
            load_end = load['end']
            length = nodes[i+1] - nodes[i]  # Length of the current element
            
            # Determine if load affects the current element
            if load_start < nodes[i+1] and load_end > nodes[i]:
                # Determine 'a' and 'b' within the current element
                if load_start < nodes[i] and load_end == nodes[i+1]:  # Case 2.
                    a1 = nodes[i]
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start < nodes[i] and load_end > nodes[i+1]:  # Case 1.
                    a1 = nodes[i]
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start < nodes[i] and load_end < nodes[i+1]:  # Case 3.
                    a1 = nodes[i]
                    b1 = load_end
                    a=0
                    b=b1-a1
                elif load_start == nodes[i] and load_end == nodes[i+1]:  # Case 4.
                    a1 = nodes[i]
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start == nodes[i] and load_end < nodes[i+1]:  # Case 5.
                    a1 = nodes[i]
                    b1 = load_end
                    a=0
                    b=b1-a1
                elif load_start == nodes[i] and load_end > nodes[i+1]:  # Case 9
                    a1 = nodes[i]
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start > nodes[i] and load_end == nodes[i+1]:  # Case 6.
                    a1 = load_start
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start > nodes[i] and load_end > nodes[i+1]:  # Case 7.
                    a1 = load_start
                    b1 = nodes[i+1]
                    a=0
                    b=b1-a1
                elif load_start > nodes[i] and load_end < nodes[i+1]:  # Case 8.
                    a1 = load_start
                    b1 = load_end
                    a=0
                    b=b1-a1
                # Calculate fixed-end moments and reactions for the segment
                UDL_moment_A = abs((w / length**2) * ((length**2 / 2) * (b**2 - a**2) - (2 * length / 3) * (b**3 - a**3) + (1 / 4) * (b**4 - a**4)))
                UDL_moment_B = abs((w / length**2) * ((length / 3) * (b**3 - a**3) - (1 / 4) * (b**4 - a**4)))
                UDL_reaction_A = abs((w * (b - a) / length) * ((length - (b + a) / 2)))
                UDL_reaction_B = abs((w * (b - a) / length) * ((b + a) / 2))

                # Adjust for load direction
                if load['direction'] == 'down':
                    UDL_reaction_A = -UDL_reaction_A
                    UDL_reaction_B = -UDL_reaction_B
                    UDL_moment_A = -UDL_moment_A
                if load['direction'] == 'up':
                    UDL_moment_B = -UDL_moment_B

                # Apply load effects to element force vector
                element_force_vector[2*i] += UDL_reaction_A
                element_force_vector[2*i+1] += UDL_moment_A
                element_force_vector[2*(i+1)] += UDL_reaction_B
                element_force_vector[2*(i+1)+1] += UDL_moment_B
        
        elif load['type'] == 'moment':
            if nodes[i] <= load['position'] < nodes[i+1]:
                # Calculate distances to the left and right nodes
                distance_to_left = load['position'] - nodes[i]
                distance_to_right = nodes[i+1] - load['position']
                total_distance = nodes[i+1] - nodes[i]

                # Distribute moment based on distances
                M_left = abs((distance_to_right / total_distance) * load['value'])
                M_right = abs((distance_to_left / total_distance) * load['value'])

                if load['direction']=='clockwise':
                    M_left=-M_left
                    M_right=-M_right
                
                element_force_vector[2*i+1] += M_left
                element_force_vector[2*(i+1)+1] += M_right

    element_force_vectors.append(element_force_vector)
    print(f"Force Vector for Element {i+1}:")
    print(np.array2string(np.round(element_force_vector, 3), separator=', '))

# Sum up individual force vectors to get the overall force vector
overall_force_vector = np.sum(element_force_vectors, axis=0)

print("Overall Force Vector [F]:")
print(np.array2string(np.round(overall_force_vector, 3), separator=', '))

print("Displacement Vector [u]:")
print(displacement_vector)

# Step 6: Update overall force vector based on support types
for i, support in enumerate(supports):
    if support == 'fixed':
        # Fixed support has reaction force and moment
        force_vector[2*i] = f'R{i+1}'  # Symbolic reaction force, do not perform arithmetic
        force_vector[2*i+1] = f'M{i+1}'  # Symbolic moment, do not perform arithmetic
    elif support == 'roller':
        # Roller support has only reaction force
        force_vector[2*i] = f'R{i+1}'  # Symbolic reaction force, do not perform arithmetic
    elif support == 'hinge':
        # Hinge support has only reaction force, no moment
        force_vector[2*i] = f'R{i+1}'  # Symbolic reaction force, do not perform arithmetic
final_force_vector=[]
# Print final force vector considering supports
print("Final Force Vector [F] after considering supports:")
for i in range(global_matrix_size):
    force_component = f"{overall_force_vector[i]} {' ' if isinstance(force_vector[i], (int, float)) else'+'+force_vector[i]}"
    final_force_vector.append(force_component)
print(final_force_vector)
global_stiffness_matrix = np.array(global_stiffness_matrix)
# Convert displacement_vector and final_force_vector to symbolic expressions
U = sp.Matrix([sp.sympify(u) if isinstance(u, str) else u for u in displacement_vector])

F = sp.Matrix([sp.sympify(f) if isinstance(f, str) else f for f in final_force_vector])

# Global stiffness matrix as a sympy Matrix
K = sp.Matrix(global_stiffness_matrix)

# The equation is F = K * U
equations = F - K * U

# Solve the system of equations
unknowns = sp.solve(equations, dict=True)

# Output the result
print("Unknown Values: ")
print(unknowns)

unknowns_obj=unknowns[0]
flat_unknowns = {str(key): value for key, value in unknowns_obj.items()}

# Output the new dictionary with string keys
print(flat_unknowns)


def process_force_vector(final_force_vector, flat_unknowns, global_matrix_size):
    def evaluate_expression(expr, unknowns):
        """Evaluate a string expression using sympy"""
        if isinstance(expr, (int, float)):
            return float(expr)
        
        expr = sp.sympify(expr)
        symbols = expr.free_symbols
        
        # If the expression is just a single symbol, look it up in unknowns
        if len(symbols) == 1 and str(list(symbols)[0]) in unknowns:
            return float(unknowns[str(list(symbols)[0])])
        
        # Otherwise, substitute known values and evaluate
        return float(expr.evalf(subs={str(s): unknowns.get(str(s), s) for s in symbols}))

    r_values = []
    m_values = []

    for i in range(0, global_matrix_size, 2):
        # Process reaction forces (R)
        r_expr = final_force_vector[i]
        r_value = evaluate_expression(r_expr, flat_unknowns)
        r_values.append(r_value)

        # Process moments (M)
        m_expr = final_force_vector[i+1]
        m_value = evaluate_expression(m_expr, flat_unknowns)
        m_values.append(m_value)

    return r_values, m_values

# Process the force vector
r_values, m_values = process_force_vector(final_force_vector, flat_unknowns, global_matrix_size)

# Process displacement vector
v_values = []
theta_values = []

for i in range(0, global_matrix_size, 2):
    v_value = displacement_vector[i]
    v_values.append(float(flat_unknowns.get(str(v_value), v_value)))

    theta_value = displacement_vector[i+1]
    theta_values.append(float(flat_unknowns.get(str(theta_value), theta_value)))

# Output the separated arrays
print("Vertical Displacements (v):", v_values)
print("Rotations (theta):", theta_values)
print("Reactions (r):", r_values)
# print("Moments (m):", m_values)

# Initialize shear_force array
shear_force = []

# Iterate over r_values to compute cumulative sum/difference for shear_force
cumulative_sum = 0
for value in r_values:
    cumulative_sum += value
    shear_force.append(cumulative_sum)

x_values = [i * element_length for i in range(num_elements + 1)]
# Plotting Vertical Displacements (v)
plt.figure(figsize=(10, 6))
plt.plot(x_values, v_values, marker='o', label="Vertical Displacement (v)")
plt.fill_between(x_values, v_values, color='skyblue', alpha=0.6)  # Shading
plt.xlabel('Length along the beam (m)')
plt.ylabel('Vertical Displacement (v)')
plt.title('Vertical Displacements along the Beam')
plt.grid(True)
plt.legend()
plt.show()

# Plotting Rotations (theta)
plt.figure(figsize=(10, 6))
plt.plot(x_values, theta_values, marker='o', label="Rotations (theta)")
plt.fill_between(x_values, theta_values, color='lightcoral', alpha=0.6)  # Shading
plt.xlabel('Length along the beam (m)')
plt.ylabel('Rotation (theta)')
plt.title('Rotations along the Beam')
plt.grid(True)
plt.legend()
plt.show()

# Plotting Reactions (r)
plt.figure(figsize=(10, 6))
plt.plot(x_values, shear_force, marker='o', label="Reactions (r)")
plt.fill_between(x_values, shear_force, color='lightgreen', alpha=0.6)  # Shading
plt.xlabel('Length along the beam (m)')
plt.ylabel('Reaction (r)')
plt.title('Reactions along the Beam')
plt.grid(True)
plt.legend()
plt.show()

# Plotting Moments (m)
# plt.figure(figsize=(10, 6))
# plt.plot(x_values, m_values, marker='o', label="Moments (m)")
# plt.fill_between(x_values, m_values, color='lightyellow', alpha=0.6)  # Shading
# plt.xlabel('Length along the beam (m)')
# plt.ylabel('Moment (m)')
# plt.title('Moments along the Beam')
# plt.grid(True)
# plt.legend()
# plt.show()