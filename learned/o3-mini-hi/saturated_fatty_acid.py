"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: Any fatty acid containing no carbon to carbon multiple bonds.
They are characterized by a carboxylic acid group and a fully saturated (i.e. no C=C or C≡C) hydrocarbon chain.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    A saturated fatty acid must contain a carboxylic acid group (-C(=O)OH) and should not have any carbon-carbon
    double or triple bonds (unsaturation) outside of the carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a carboxylic acid group.
    # SMARTS pattern: C(=O)[O;H,-] matches a carbon with a double-bonded oxygen and an oxygen with either a hydrogen or a negative charge.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for any carbon-carbon multiple bonds (double or triple).
    # We iterate over bonds in the molecule; allow the carbonyl C=O but any other C=C or C≡C is disallowed.
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            # Get the two atoms involved in the bond.
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # If both atoms are carbon, then check if this bond is part of the carboxyl group.
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                # Exclude the carbonyl in the acid group by checking if one carbon is connected to an oxygen in a C=O.
                # Get neighbors of the first carbon; if one neighbor is an oxygen with a double bond then assume it's the carbonyl.
                # Since both carbons are from a C=C in many cases, it's safer to simply return unsaturated if any C=C (outside of the acid group) exists.
                return False, "Contains carbon–carbon multiple bond(s); not a saturated fatty acid"
    
    return True, "Saturated fatty acid: contains carboxylic acid group and no carbon–carbon multiple bonds"