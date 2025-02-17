"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is defined as a compound whose main (largest)
            organic fragment contains at least one carbon–chlorine bond.
            
Improvements in this version:
  - Split the molecule into fragments and work on the largest one.
  - Ensure that the largest fragment contains only allowed elements.
  - For sufficiently large fragments (≥10 heavy atoms) require a modest organic character:
       • carbon fraction must be at least 0.35  (relaxed from 0.40)
       • heteroatom-to-carbon ratio (non-C/C) must be at most 2.5.
  - Then, instead of blindly matching a [#6]-[Cl] SMARTS pattern, we “validate” each C–Cl bond.
    For a bond to count:
       • If the fragment is very small (<5 heavy atoms) any C-Cl counts.
       • Otherwise, the chlorine (Cl) must be bound to a carbon that has at least one other
         carbon neighbor (i.e. is attached to some hydrocarbon chain).        
  - This heuristic helps filter out cases where a chlorine is attached on a peripheral position 
    in a system that is otherwise not “organic‐rich.”
    
Note:
  This is still a heuristic approach. Some borderline cases might be misclassified.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound.
    
    For our purposes the molecule is classified as an organochlorine compound if its largest
    fragment is built only from allowed ('organic') atoms, is "organic enough" (when large)
    and that fragment contains at least one validated C–Cl bond.
    
    A validated C–Cl bond is one where:
      - For a small fragment (<5 heavy atoms) any C–Cl bond counts.
      - For larger fragments the carbon atom attached to chlorine must have at least one other
        neighboring carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise.
        str: Reason for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split into fragments and take the largest (by # of heavy atoms).
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found"
    largest_frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    # Check that the largest fragment contains at least one carbon.
    if not any(atom.GetAtomicNum() == 6 for atom in largest_frag.GetAtoms()):
        return False, "Largest fragment contains no carbon atoms"
    
    # Allowed elements typical for organic compounds:
    allowed_atomic_nums = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}
    for atom in largest_frag.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non‐organic element: {atom.GetSymbol()}"
    
    # For larger fragments enforce organic character.
    heavy_atoms = [atom for atom in largest_frag.GetAtoms() if atom.GetAtomicNum() > 1]
    num_heavy = len(heavy_atoms)
    num_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    
    # If the fragment is sufficiently large, check carbon fraction and heteroatom ratio.
    if num_heavy >= 10:
        carbon_fraction = num_carbon / num_heavy
        if carbon_fraction < 0.35:
            return False, f"Organic fraction too low (carbon fraction = {carbon_fraction:.2f})"
        hetero_ratio = (num_heavy - num_carbon) / num_carbon if num_carbon > 0 else float('inf')
        if hetero_ratio > 2.5:
            return False, f"Excessive heteroatom content (non-C/C ratio = {hetero_ratio:.2f})"
    
    # Now check for at least one valid C–Cl bond.
    # Iterate over all chlorine atoms in the largest fragment.
    valid_ccl_found = False
    for atom in largest_frag.GetAtoms():
        if atom.GetAtomicNum() == 17:  # chlorine
            for nbr in atom.GetNeighbors():
                # Check if neighbor is carbon.
                if nbr.GetAtomicNum() == 6:
                    # For small fragments, accept any C-Cl.
                    if num_heavy < 5:
                        valid_ccl_found = True
                        break
                    # For larger fragments, check that this carbon has at least one other carbon neighbor.
                    carbon_neighbor_count = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom.GetIdx())
                    if carbon_neighbor_count >= 1:
                        valid_ccl_found = True
                        break
            if valid_ccl_found:
                break

    if not valid_ccl_found:
        return False, "Does not contain any validated carbon–chlorine bonds in the main organic fragment"
    
    return True, "Contains at least one validated carbon–chlorine bond in the main organic fragment"


# Example usage to test several SMILES strings:
if __name__ == "__main__":
    test_smiles = [
        "ClCCBr",  # 1-bromo-2-chloroethane
        "C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)CN3C=NN=N3)O[C@H]1CN(C)CC4=CC(=C(C=C4)Cl)Cl)[C@@H](C)CO",
        "ClC(Cl)(Cl)C(c1ccccc1)c1ccccc1",
        "CCOC(=O)Nc1ccc(SCC2COC(Cn3ccnc2)c2ccc(Cl)cc2Cl)cc1",
        "Clc1ccc(Cl)c(c1)-c1cc(Cl)ccc1Cl",
        "C1CCC(C1)CN2C[C@H](COC[C@H]3[C@H]2CC[C@@H](O3)CC(=O)NCC4=CC=C(C=C4)Cl)O",
        "Clc1c(Cl)c(Cl)c(-c2ccccc2)c(Cl)c1Cl",
        "[C@H]1(C([C@@H]2CC[C@]1(C2)C)(C)C)NC(=O)C3=NN(C(=C3)C4=CC(=C(Cl)C=C4)C)CC5=CC=C(C=C5)C",
        "OC(=O)\\C=C\\c1ccc(Cl)cc1",
        "Oc1ccc(Cl)cc1-c1cc(-c2cc(OCc3ccccc3)cc(OCc3ccccc3)c2)c(C#N)c(=O)[nH]1",
        # Several cases that were problematic in the previous attempt:
        "ClC[C@H]1N[C@H](C(=O)N[C@H](C(=O)N2[C@H](C(=O)N3[C@H](C(=O)O)C[C@@H](C3)C)CCCC2)C(C)C)CC[C@@H]1OC",  # false positive example (should be False)
        "C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O",  # false positive example (should be False)
        "C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]3[C@@H](O2)COC[C@@H](CN3C(=O)NC4=CC(=CC=C4)Cl)O",  # false positive example (should be False)
        "ClC1=C(OC)C(C(=O)C=2NC(Cl)=C(C2)Cl)=CC(=C1)Cl",  # false positive example (should be False)
        "OC(=O)c1c(Cl)c(Cl)c(Cl)c(Cl)c1-c1c2cc(Br)c(O)c(Br)c2oc2c(Br)c(=O)c(Br)cc12",  # false positive example (should be False)
        "C(C(CO)O)N1C=C(I)C(C(=C1)I)=O",  # iopydol: no C-Cl, should be False.
        "C1(=NC(=C2C(=N1)N(C=N2)[C@@H]3O[C@@H]([C@H]([C@H]3O)O)COP(OP(C(P(O)(=O)O)(Cl)Cl)(=O)O)(=O)O)NCCSC)SCCC(F)(F)F",  # cangrelor false negative: now should pass if it has a valid C-Cl
        "N1C2SC(=NN2C(=N1)C=3SN=NC3C)C(Cl)(Cl)Cl",  # trichloromethyl triazolothiadiazole (should pass if organic enough) 
        "OC(=O)C=Cc1ccccc1",  # no Cl so False.
        "CC1=NN=C(O1)N=C2C(=NSS2)Cl",  # should detect C-Cl if organic enough.
        "FC(F)(F)[C@@H](Cl)Br",  # (S)-halothane, should be True.
        "ClC(Br)Br",  # Chlorodibromomethane, should be True.
        "O=C(O)CC(Cl)C(O)=O",  # 2-chlorosuccinic acid, should be True.
    ]
    for smi in test_smiles:
        res, reason = is_organochlorine_compound(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")