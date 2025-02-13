"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:35655 lactol

A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group 
to an aldehydic or ketonic carbonyl group. They are thus 1-oxacycloalkan-2-ols or 
unsaturated analogues, typically 5- or 6-membered rings.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and get tautomers
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    tautomers = [Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True))]
    tautomers.extend(Chem.AddHs(mol).Salts().EnumerateTautomers())
    
    # Check each tautomer
    for tautomer in tautomers:
        ssr = rdchem.GetSSSR(tautomer)
        for ring in ssr:
            if len(ring) in [5, 6]:  # Ring size 5 or 6
                ring_atoms = [tautomer.GetAtomWithIdx(idx) for idx in ring]
                ring_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() != 1]  # Exclude hydrogens
                oxygens = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
                if len(oxygens) == 1:  # Single ring oxygen
                    oxygen = oxygens[0]
                    if oxygen.GetDegree() == 2:  # Cyclic ether oxygen
                        neighbors = [tautomer.GetAtomWithIdx(nbr) for nbr in oxygen.GetNeighbors()]
                        carbonyls = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 0]
                        if len(carbonyls) == 1:  # Adjacent carbonyl carbon
                            carbonyl = carbonyls[0]
                            if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in carbonyl.GetBonds()):
                                # Carbonyl carbon is part of a C=O group
                                return True, "Molecule contains a lactol structure"
                            
    return False, "No lactol structure found"