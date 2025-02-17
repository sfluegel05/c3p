"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: Polypeptide, defined as a peptide containing ten or more amino acid residues.
This improved version uses RDKit to iterate over bonds and only count those amide bonds that 
appear to be part of a peptide backbone. For each bond between a carbon and a nitrogen, we require 
that (a) the carbon has at least one neighbor with a double bond that is oxygen (the carbonyl) and 
(b) the nitrogen is linked to at least one other carbon atom. We then determine if the peptide chain 
appears “linear” (free -NH2 and -COOH termini detected) or cyclic/blocked.
Note: This is still a heuristic and may mis‐count residues in some edge‐cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (≥10 amino acid residues) based on its SMILES string.
    Uses a heuristic by counting peptide bonds that appear on a peptide backbone.
    
    Algorithm:
      – Parse the molecule.
      – For each bond, if one end is carbon and the other nitrogen, check that:
            (a) The carbon is part of a group C(=O)... (i.e. at least one doublebond to O)
            (b) The nitrogen has at least one other bond to carbon (i.e. not just an acyl NH).
         Such bonds are counted as peptide bonds.
      – Then, check for free terminal groups:
            If a free primary amine ([NH2]) and a free carboxylic acid (C(=O)[OH]) are found,
            assume “linear” so that residue count = peptide bond count + 1.
         Otherwise, assume cyclic or blocked termini so that residue count equals peptide bond count.
      – Return True if the estimated residue count ≥ 10.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a polypeptide, False otherwise.
        str: Reason for classification.
    """
    # Convert SMILES into a molecule using RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    peptide_bond_count = 0

    # Loop over bonds to count only those that look like peptide bonds.
    for bond in mol.GetBonds():
        # We want bonds between a carbon and a nitrogen (order can be C-N or N-C).
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify which atom is C and which is N
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7:
            carbon_atom = a1
            nitrogen_atom = a2
        elif a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6:
            nitrogen_atom = a1
            carbon_atom = a2
        else:
            continue  # Skip bonds not between C and N

        # Check that the carbon_atom has at least one neighbor oxygen involved in a double bond.
        carbon_has_carbonyl = False
        for nbr in carbon_atom.GetNeighbors():
            # Make sure the neighbor is oxygen and the connecting bond is a double bond.
            if nbr.GetAtomicNum() == 8:
                b = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                if b is not None and b.GetBondTypeAsDouble() >= 2.0:
                    carbon_has_carbonyl = True
                    break
        if not carbon_has_carbonyl:
            continue

        # Check that the nitrogen_atom is connected to at least one other carbon (besides the found carbon)
        nitrogen_has_c_neighbor = False
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetIdx() == carbon_atom.GetIdx():
                continue  # skip our carbon partner
            if nbr.GetAtomicNum() == 6:
                nitrogen_has_c_neighbor = True
                break
        if not nitrogen_has_c_neighbor:
            continue

        # If both conditions met, count this bond as a peptide bond.
        peptide_bond_count += 1

    # Now estimate if the peptide looks linear by searching for free terminal groups.
    # For a free N-terminus use a primary amine, and for a free C-terminus use a carboxylic acid.
    free_amine_pattern = Chem.MolFromSmarts("[NH2]")
    free_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    has_free_amine = mol.HasSubstructMatch(free_amine_pattern)
    has_free_acid = mol.HasSubstructMatch(free_acid_pattern)
    
    if has_free_amine and has_free_acid:
        # For linear peptides, residue count = number of peptide bonds + 1.
        residue_count = peptide_bond_count + 1
        peptide_type = "linear"
    else:
        # For cyclic or blocked peptides, residue count ~ number of peptide bonds.
        residue_count = peptide_bond_count
        peptide_type = "cyclic or with blocked termini"
    
    if residue_count >= 10:
        return True, f"Detected {residue_count} amino acid residues (peptide appears {peptide_type})."
    else:
        return False, f"Detected only {residue_count} amino acid residues (need at least 10)."

# Example usage:
if __name__ == "__main__":
    # Test examples – a few examples are provided below.
    examples = {
        "nisin": "CC[C@H](C)[C@H](N)C(=O)N\\C(=C/C)C(=O)N[C@@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(=C)NC(=O)[C@@H](NC1=O)[C@@H](C)CC)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@@H]2CCCN2C1=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]1[C@H](C)SC[C@H](NC(=O)CNC(=O)[C@H](CCSC)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](C)NC(=O)CNC1=O)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(N)=O",
        "Ac-Asp-N(6)-[2-(trifluoromethyl)benzoyl]-KATIGFEVQEE": "CC[C@H](C)[C@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCNC(=O)c1ccccc1C(F)(F)F)NC(=O)[C@H](CC(O)=O)NC(C)=O)[C@@H](C)O)C(=O)NCC(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CCC(O)=O)C(O)=O"
    }
    
    for name, sm in examples.items():
        status, reason = is_polypeptide(sm)
        print("NAME:", name)
        print("Is polypeptide?", status)
        print("Reason:", reason, "\n")