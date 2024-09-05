import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import re

# Set the print options for NumPy arrays to show 3 decimal places
# np.set_printoptions(precision=3, suppress=True)

# Step 1: Get user input
beam_length = float(input("Enter the total length of the beam (in meters): "))
num_elements = int(input("Enter the number of elements you want: "))
num_nodes = num_elements + 1

# Determine node positions
element_length = beam_length / num_elements
nodes = np.linspace(0, beam_length, num_nodes)

supports = ["none"] * num_nodes

# Get support types
intermediate_support = input(f"Are there any intermediate support(s)? (Yes/No): ").strip().lower()

# Get start and end supports
supports[0] = input(f"Enter the type of support at the start (node 1) (Fixed, Roller, Hinge, None): ").strip().lower()
supports[-1] = input(f"Enter the type of support at the end (node {num_nodes}) (Fixed, Roller, Hinge, None): ").strip().lower()

if intermediate_support == "yes":
    num_supports = int(input("Enter the number of intermediate supports: "))
    
    for _ in range(num_supports):
        support_distance = float(input("Enter the distance of the intermediate support from the start (in meters): "))
        support_type = input("Enter the type of support (Fixed, Roller, Hinge): ").strip().lower()

        # Find the closest node
        closest_node = np.argmin(np.abs(nodes - support_distance))

        if nodes[closest_node] > support_distance and closest_node > 0:
            closest_node -= 1
        
        supports[closest_node] = support_type

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
print(f"Nodes (equal division): {nodes}")

# Step 2: Calculate element lengths
element_lengths = np.diff(nodes)
# print(f"Element Lengths: {element_lengths}")

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
    # print(f"Stiffness Matrix for Element {i+1}:")
    # print(np.array2string(k, separator=', '))  # Improved matrix printing with 3 decimals and converted to string for separation by comma

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
print(np.array2string(global_stiffness_matrix, separator=', '))  

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
                a = load['position'] - nodes[i]
                b = nodes[i+1] - load['position']
                length = nodes[i+1] - nodes[i]
                
                # Updated Reaction Formulas
                F_left = abs(load['value'] * b**2 / length**3 * (length + 2*a))
                F_right = abs(load['value'] * a**2 / length**3 * (length + 2*b))
                
                # Updated Moment Formulas
                M_left = abs(load['value'] * a * b**2 / length**2)
                M_right = abs(load['value'] * a**2 * b / length**2)
                
                if load['direction'] == 'down':
                    F_left = -F_left
                    F_right = -F_right
                    M_left = -M_left
                elif load['direction'] == 'up':
                    M_right = -M_right

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
                    a=load_start-nodes[i]
                    b=a+(b1-a1)
                elif load_start > nodes[i] and load_end > nodes[i+1]:  # Case 7.
                    a1 = load_start
                    b1 = nodes[i+1]
                    a=load_start-nodes[i]
                    b=a+(b1-a1)
                elif load_start > nodes[i] and load_end < nodes[i+1]:  # Case 8.
                    a1 = load_start
                    b1 = load_end
                    a=load_start-nodes[i]
                    b=a+(b1-a1)
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
    # print(f"Force Vector for Element {i+1}:")
    # print(np.array2string(element_force_vector, separator=', '))

# Sum up individual force vectors to get the overall force vector
overall_force_vector = np.sum(element_force_vectors, axis=0)

# print("Overall Force Vector [F]:")
# print(np.array2string(overall_force_vector, separator=', '))

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
print("\nUnknown Values: \n")
# print(unknowns)

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
# print("Vertical Displacements (v):", v_values)
# print("Rotations (theta):", theta_values)
# print("Reactions (r):", r_values)
# print("Moments (m):", m_values)

# Initialize arrays for left and right shear forces
sf_left = [0] * len(nodes)
sf_right = [0] * len(nodes)

def calculate_load_at_point(x):
    total_load = 0
    print(f"Calculating load at x = {x}:")
    for load in loads:
        if load['type'] == 'point' and load['position'] < x:
            load_value = -load['value'] if load['direction'] == 'down' else load['value']
            print(f"  Point load at position {load['position']} contributes: {load_value}")
            total_load += load_value
        elif load['type'] == 'udl':
            start = min(x, max(load['start'], 0))
            end = min(x, load['end'])
            if end > start:
                length = end - start
                udl_value = load['value'] * length * (-1 if load['direction'] == 'down' else 1)
                print(f"  UDL from {start} to {end} contributes: {udl_value} (Length = {length})")
                total_load += udl_value
    print(f"  Total load at x = {x}: {total_load}")
    return total_load

# Calculate shear forces
cumulative_load = 0
cumulative_reaction = 0

for i, node in enumerate(nodes):
    print(f"\nNode {i+1} at position {node}:")
    
    # Calculate cumulative load up to this point for sf_left
    cumulative_load = calculate_load_at_point(node)
    
    # Calculate sf_left
    sf_left[i] = cumulative_reaction + cumulative_load
    print(f"  Cumulative load up to node {node}: {cumulative_load}")
    print(f"  Cumulative reaction up to node {node}: {cumulative_reaction}")
    print(f"  SF Left before adding support: {sf_left[i]}")

    # Add support reaction to cumulative reaction if present at this node
    if supports[i] in ['fixed', 'roller', 'hinge']:
        print(f"  Adding support reaction {r_values[i]} at node {node}")
        cumulative_reaction += r_values[i]
    
    # Calculate sf_right
    sf_right[i] = cumulative_reaction + cumulative_load
    print(f"  SF Right before checking point loads: {sf_right[i]}")
    
    # Check for point load exactly at this node and add to sf_right
    for load in loads:
        if load['type'] == 'point' and load['position'] == node:
            load_value = -load['value'] if load['direction'] == 'down' else load['value']
            print(f"  Point load exactly at node {node} contributes: {load_value}")
            sf_right[i] += load_value
    
    print(f"  Updated SF Right: {sf_right[i]}")

# Print sf_right last value
print(f"\nSF_right last value: {sf_right[-1]}")

# Plotting Shear Force
plt.figure(figsize=(10, 6))

# Plot shear forces with correct connections
for i in range(len(nodes) - 1):
    plt.plot([nodes[i], nodes[i]], [sf_left[i], sf_right[i]], 'b-')  # Vertical line at each node
    plt.plot([nodes[i], nodes[i+1]], [sf_right[i], sf_left[i+1]], 'r-')  # Connecting line between nodes

# Plot markers for left and right shear forces
plt.plot(nodes, sf_left, 'bo', label="SF Left", markersize=4)
plt.plot(nodes, sf_right, 'ro', label="SF Right", markersize=4)

# Fill between shear force lines
for i in range(len(nodes) - 1):
    plt.fill_between([nodes[i], nodes[i+1]], 
                     [sf_right[i], sf_left[i+1]], 
                     color='lightblue', alpha=0.8)

plt.xlabel('Length along the beam (m)')
plt.ylabel('Shear Force (kN)')
plt.title('Shear Force Diagram')
plt.grid(True)
plt.legend()
plt.show()


cumulative_bm_left = []
cumulative_bm_right = []

# Function to find the closest node
def find_closest_node(position, nodes):
    return np.argmin(np.abs(nodes - position))

# Calculate bending moments
for i, node in enumerate(nodes):
    print(f"\nNode {i + 1} at position {node}:")
    
    # Initialize arrays for bending moments at this node
    bm_left = 0
    bm_right = 0
    
    print(f"Initial BM Left: {bm_left}, BM Right: {bm_right}")
    
    # Add moments due to external moments
    for load in loads:
        if load['type'] == 'moment':
            closest_node_idx = find_closest_node(load['position'], nodes)
            if closest_node_idx <(i):  # Only consider loads to the left
                if load['direction'] == 'clockwise':
                    print(f"Moment load at node {closest_node_idx + 1}: Adding {load['value']} clockwise to BM Right and BM Left")
                    bm_left+=abs(load['value'])
                    bm_right += abs(load['value'])
                else:
                    print(f"Moment load at node {closest_node_idx + 1}: Subtracting {load['value']} counterclockwise from BM Right and BM Left")
                    bm_left-=abs(load['value'])
                    bm_right -= abs(load['value'])
                print(f"Updated BM Right: {bm_right}")
            elif closest_node_idx == i:
                if load['direction'] == 'clockwise':
                    print(f"Moment load at node {closest_node_idx + 1}: Adding {load['value']} clockwise to BM Right and BM Left")
                    bm_right += abs(load['value'])
                else:
                    print(f"Moment load at node {closest_node_idx + 1}: Subtracting {load['value']} counterclockwise from BM Right")
                    bm_right -= abs(load['value'])

    
    # Add moments due to reaction forces from supports to the left
    for j in range(i):  # Only consider supports to the left of or at the current node
        support = supports[j]
        if support in ['fixed', 'roller', 'hinge']:
            distance = abs(nodes[j] - node)
            reaction_force = r_values[j]
            print(f"Support {j + 1}: Reaction force = {reaction_force}, Distance = {distance}")
            moment = reaction_force * distance
            print(f"Moment due to support = {reaction_force} * {distance} = {moment}")
            bm_left += moment
            bm_right += moment
            print(f"Updated BM Left: {bm_left}, BM Right: {bm_right}")

    # Add moments due to point loads to the left
    for load in loads:
        if load['type'] == 'point' and load['position'] < node:
            distance = node - load['position']
            moment = -load['value'] * distance if load['direction'] == 'down' else load['value'] * distance
            print(f"Point load: Distance = {distance}, Moment = {load['value']} * {distance} = {moment}")
            bm_left += moment
            bm_right += moment
            print(f"Updated BM Left: {bm_left}, BM Right: {bm_right}")

    # Add moments due to UDLs to the left
    for load in loads:
        if load['type'] == 'udl' and load['start'] < node:
            if node > load['end']:  # Node is beyond UDL
                length = load['end'] - load['start']
                distance = node - load['end']
            else:  # Node is within the UDL
                length = node - load['start']
                distance = 0
            udl_moment = load['value'] * length * ((length / 2) + distance)
            direction_multiplier = -1 if load['direction'] == 'down' else 1
            udl_moment *= direction_multiplier
            print(f"UDL: Length = {length}, Distance = {distance}, UDL Moment = {udl_moment}")
            bm_left += udl_moment
            bm_right += udl_moment
            print(f"Updated BM Left: {bm_left}, BM Right: {bm_right}")

    # Add induced moments at fixed supports at this node
    
    for k in range(i + 1):
        if supports[k] == 'fixed':
            if k < i:
                print(f"Fixed support at node {k + 1}: Adding induced moment {m_values[k]} to both BM Left and BM Right")
                bm_left -= m_values[k]
                bm_right -= m_values[k]
            else:  # k == i
                print(f"Fixed support at node {i + 1}: Adding induced moment {m_values[k]} to BM Right only")
                bm_right -= m_values[k]
            print(f"Updated BM Left: {bm_left}, BM Right: {bm_right}")

    # Append bending moments to cumulative arrays
    cumulative_bm_left.append(bm_left)
    cumulative_bm_right.append(bm_right)
    print(f"Final BM Left: {bm_left}, BM Right: {bm_right}")

# Print bending moments for all nodes
print("Bending Moments (Left):", cumulative_bm_left)
print("Bending Moments (Right):", cumulative_bm_right)

# Plotting Bending Moment Diagram
plt.figure(figsize=(10, 6))

# Plot bending moments with correct connections
for i in range(len(nodes) - 1):
    plt.plot([nodes[i], nodes[i]], [cumulative_bm_left[i], cumulative_bm_right[i]], 'b-')  # Vertical line at each node
    plt.plot([nodes[i], nodes[i+1]], [cumulative_bm_right[i], cumulative_bm_left[i+1]], 'r-')  # Connecting line between nodes

# Plot markers for left and right bending moments
plt.plot(nodes, cumulative_bm_left, 'bo', label="BM Left", markersize=4)
plt.plot(nodes, cumulative_bm_right, 'ro', label="BM Right", markersize=4)

# Fill between bending moment lines
for i in range(len(nodes) - 1):
    plt.fill_between([nodes[i], nodes[i+1]], 
                     [cumulative_bm_right[i], cumulative_bm_left[i+1]], 
                     color='lightblue', alpha=0.8)

plt.xlabel('Length along the beam (m)')
plt.ylabel('Bending Moment (kNm)')
plt.title('Bending Moment Diagram')
plt.grid(True)
plt.legend()
plt.show()

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

