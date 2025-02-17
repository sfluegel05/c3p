"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: Pterocarpans – Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton (and substituted derivatives).
They are 3,4-dihydro derivatives of coumestans and generally resemble isoflavanoids, often with antibiotic activity.
"""

from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    
    The classification is based on identification of a pterocarpan core scaffold.
    Pterocarpans are defined here as molecules that show the benzofurochromene core. Because 
    drawings may vary (for example, in many pterocarpans the coumestan core is partially reduced
    into a 3,4-dihydro species), we use two SMARTS patterns:
    
      1. A “fully aromatic” benzofuran pattern: "c1ccc2occc2c1"
         (this may catch some coumestan-like structures)
         
      2. A pattern that reflects the dihydro (saturated) variant in which two ring carbons are sp3: 
         "O1[CX4][CX4]Oc2ccccc2C1"
         
    If either pattern is found, we assume that the molecule contains a pterocarpan-like scaffold.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule contains a pterocarpan-like core, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string into an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns:
    # Pattern 1: fully aromatic benzofuran (relaxed coumestan core)
    arom_benzofuran = Chem.MolFromSmarts("c1ccc2occc2c1")
    if arom_benzofuran is None:
        return False, "Internal error: could not parse aromatic benzofuran SMARTS"
    
    # Pattern 2: dihydrobenzofuran core (with two saturated carbons, using CX4)
    dihydro_pattern = Chem.MolFromSmarts("O1[CX4][CX4]Oc2ccccc2C1")
    if dihydro_pattern is None:
        return False, "Internal error: could not parse dihydro benzofuran SMARTS"
    
    # Check for either pattern being present in the molecule.
    if mol.HasSubstructMatch(arom_benzofuran):
        return True, "Molecule contains an aromatic benzofuran core scaffold (pterocarpan-like)"
    
    if mol.HasSubstructMatch(dihydro_pattern):
        return True, "Molecule contains a dihydrobenzofuran core scaffold (pterocarpan-like)"
    
    # If no pattern matched, it is assumed not to be a pterocarpan.
    return False, "Pterocarpan core (benzofurochromene skeleton) not found"

    
# Example usage (for debugging/testing purposes):
if __name__ == "__main__":
    # List several example SMILES strings for pterocarpan compounds.
    test_smiles = [
        "O1C2C(C3=C1C=C(O)C=C3)COC4=C2C=CC(O)=C4O",                   # 3,4,9-Trihydroxypterocarpan
        "O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4",           # 2,8-Dihydroxy-3,9,10-trimethoxypterocarpan
        "COc1c(CC=C(C)C)c(O)cc2OC[C@@H]3[C@@H](Oc4cc(O)ccc34)c12"         # edudiol
    ]
    
    for smi in test_smiles:
        res, reason = is_pterocarpans(smi)
        print(f"SMILES: {smi}\nResult: {res}, Reason: {reason}\n")