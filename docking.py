import os

# Parameters
receptor = "DNMT1_H_cleaned.pdbqt"
ligand_dir = "/home/hp/Desktop/NIsha/ST.Xavier/Final_Project/GSE248260/Docking/ligands"
output_dir = "/home/hp/Desktop/NIsha/ST.Xavier/Final_Project/GSE248260/Docking/results"
center = (-11.023, 2.851, 29.938)
size = (75.222, 75.222, 75.222)


# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Loop through ligand files
for i in range(1, 5):
    ligand_file = f"ligand{i}.pdbqt"
    ligand_path = os.path.join(ligand_dir, ligand_file)
    out_path = os.path.join(output_dir, f"ligand{i}_out.pdbqt")
    log_path = os.path.join(output_dir, f"ligand{i}.log")
    
    cmd = (
        f"vina --receptor {receptor} "
        f"--ligand {ligand_path} "
        f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} "
        f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} "
        f"--out {out_path} > {log_path} 2>&1"
    )


    print(f"Running docking for: ligand{i}")
    os.system(cmd)
