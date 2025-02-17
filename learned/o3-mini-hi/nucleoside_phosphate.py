"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
#!/usr/bin/env python
"""
Classifies: Nucleoside phosphate
A nucleobase-containing molecular entity that is a nucleoside in which 
one or more of the sugar hydroxy groups has been converted into a mono- or poly‐phosphate.
The term includes both nucleotides and non‐nucleotide nucleoside phosphates.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    
    Strategy:
    1. Parse the SMILES using RDKit.
    2. Check that the molecule contains at least one phosphorus atom.
    3. Try to locate a phosphate moiety by matching a typical substructure OP(=O)(O) 
       (this catch many cases where a hydroxyl group on a sugar is phosphorylated).
    4. Identify a nucleobase by looking for at least one aromatic ring that contains 
       two or more nitrogen atoms.
    5. Identify a sugar moiety by scanning the ring systems for a five-membered ring 
       with one oxygen and four carbons (a proxy for furanose).
       
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a nucleoside phosphate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for presence of at least one phosphorus atom
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atoms found; unlikely to be a nucleoside phosphate"
    
    # 2. Search for phosphate substructure.
    # Many examples use the fragment "OP(=O)(O)" for the phosphodiester/phosphate link.
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate linkage (OP(=O)(O)) found in the molecule"
    
    # 3. Identify nucleobase by scanning aromatic rings that contain >=2 nitrogen atoms.
    ring_info = mol.GetRingInfo()
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        # Consider only rings that are 5- or 6-membered as common in nucleobases.
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if all atoms in ring are aromatic
        if not all(atom.GetIsAromatic() for atom in atoms):
            continue
        # Count nitrogen atoms
        n_nitrogen = sum(1 for atom in atoms if atom.GetAtomicNum() == 7)
        if n_nitrogen >= 2:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No aromatic nucleobase ring (with ≥2 nitrogen atoms) found"
    
    # 4. Identify a sugar ring 
    # As a proxy we look for a 5-membered ring that contains exactly one oxygen (the ring heteroatom) 
    # and four carbons. (Ribose or deoxyribose rings generally follow this composition.)
    sugar_found = False
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_oxygen = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
            num_carbon = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
            # A typical furanose has 1 oxygen and 4 carbons.
            if num_oxygen == 1 and num_carbon == 4:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No suitable sugar ring (5-membered with 1 oxygen and 4 carbons) found"
    
    # All tests passed; we classify the molecule as a nucleoside phosphate.
    return True, "Molecule contains a phosphate group, nucleobase and sugar ring indicative of nucleoside phosphate"

# Example usage:
if __name__ == "__main__":
    # Example: 2'-deoxy-5-methyl-5'-cytidylic acid
    test_smiles = "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N"
    result, reason = is_nucleoside_phosphate(test_smiles)
    print(result, reason)