#----------------------------------------#

# Importing for script
import re
import math

#----------------------------------------#

### Functions for reading and writing parameters ###

# Read number of cells from logs/log.checkMesh file
def read_num_cells(log_path="./log.checkMesh"):
    with open(log_path, "r") as f:
        inside_mesh_stats = False
        for line in f:
            if "Mesh stats" in line:
                inside_mesh_stats = True
            elif inside_mesh_stats:
                match = re.match(r"\s*cells:\s+(\d+)", line)
                if match:
                    return int(match.group(1))
                # End of mesh stats block
                if line.strip() == "":
                    break
    raise ValueError("Could not find 'cells:' line under 'Mesh stats' in log file.")

# Write updated nCellsInCoarsestLevel to fvSolution
def update_fvsolution(file_path, new_value):
    updated_lines = []
    pattern = re.compile(r"\s*nCellsInCoarsestLevel\s+\d+;")

    with open(file_path, "r") as f:
        for line in f:
            if pattern.match(line):
                updated_lines.append(f"            nCellsInCoarsestLevel {new_value};\n")
            else:
                updated_lines.append(line)

    with open(file_path, "w") as f:
        f.writelines(updated_lines)

# Main logic
if __name__ == "__main__":
   
    log_file = "./log.checkMesh"
    fvsolution_path = "system/fvSolution"

    num_cells = read_num_cells(log_file)
    coarse_cells = round(math.sqrt(num_cells))
    
    print(f"Detected {num_cells} cells. Setting nCellsInCoarsestLevel to {coarse_cells}.")
    
    update_fvsolution(fvsolution_path, coarse_cells)
