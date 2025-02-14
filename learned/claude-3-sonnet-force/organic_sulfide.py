"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:33877 organic sulfide
Definition: Compounds having the structure RSR (R =/= H). Such compounds were once called thioethers.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, FragmentMatcher

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    
    # Check each sulfur atom
    for sulfur in sulfur_atoms:
        # Exclude thiols (-SH) and thiocyanates (-SCN)
        if any(atom.GetSmarts() == "[SH]" for atom in sulfur.GetNeighbors()) or sulfur.IsInRingOfSize(3):
            continue
        
        # Check for RSR or R=SR structure
        neighbors = [nbr for nbr in sulfur.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(neighbors) == 2:
            if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in sulfur.GetBonds()):
                return True, "Contains sulfur atom with double bond and non-hydrogen substituent"
            else:
                return True, "Contains sulfur atom with two non-hydrogen substituents"
        
        # Check for cyclic sulfides
        if len(neighbors) == 1:
            ring_atoms = mol.GetRingInfo().AtomRings()
            for ring in ring_atoms:
                if sulfur.GetIdx() in ring:
                    ring_bonds = [mol.GetBondBetweenAtoms(ring[i], ring[i-1]).GetBondType() for i in range(len(ring))]
                    if Chem.BondType.DOUBLE in ring_bonds:
                        return True, "Contains sulfur atom in a cyclic sulfide structure"
    
    # Check for aromatic sulfides
    fmatcher = FragmentMatcher.FragmentMatcher()
    if fmatcher.AromaticRings(mol).count('s') > 0:
        return True, "Contains aromatic sulfide ring system"
    
    return False, "Does not match the structure RSR (R =/= H)"

# Example usages
print(is_organic_sulfide("CNC(=O)ON=C(C)SC"))  # True, methomyl
print(is_organic_sulfide("S(C[C@@H](C(=O)O)N)C[C@H](C(O)=O)N"))  # True, meso-lanthionine
print(is_organic_sulfide("CSC1=N[C@](C)(C(=O)N1Nc1ccccc1)c1ccccc1"))  # True, fenamidone
print(is_organic_sulfide("Cn1cnc(c1Sc1ncnc2nc[nH]c12)[N+]([O-])=O"))  # True, azathioprine
print(is_organic_sulfide("SC(S)=N"))  # False, carbonimidodithioic acid (sulfur atom has only one substituent)
print(is_organic_sulfide("N[C@@H](CSCC(O)=O)C(O)=O"))  # False, S-carboxymethyl-L-cysteine (thiol group)