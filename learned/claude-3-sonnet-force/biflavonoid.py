"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: CHEBI:51811 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for various flavonoid monomers
    flavone_pattern = Chem.MolFromSmarts("c1cc(=O)c2c(O)cccc2o1")
    flavanone_pattern = Chem.MolFromSmarts("C1Cc2c(ccc(O)c2O)O1")
    flavan3ol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)c2C(=O)Cc3ccccc3O2")
    
    # Identify potential flavonoid monomers in the molecule
    monomers = []
    for pattern in [flavone_pattern, flavanone_pattern, flavan3ol_pattern]:
        matches = mol.GetSubstructMatches(pattern)
        monomers.extend([list(match) for match in matches])
    
    # Check if at least two distinct flavonoid monomers are present
    if len(set([tuple(sorted(monomer)) for monomer in monomers])) < 2:
        return False, "Fewer than two distinct flavonoid monomers found"
    
    # Check for linker groups between monomers
    linker_atoms = set()
    for i, monomer1 in enumerate(monomers):
        for j, monomer2 in enumerate(monomers):
            if i != j:
                for atom1 in monomer1:
                    for atom2 in monomer2:
                        bond = mol.GetBondBetweenAtoms(atom1, atom2)
                        if bond:
                            linker_atoms.add(atom1)
                            linker_atoms.add(atom2)
    
    # Check if linker atoms are consistent with biflavonoid linkers
    allowed_linkers = [Chem.MolFromSmarts("[#6]-[#6]"), Chem.MolFromSmarts("[#6]-[#8]-[#6]")]
    for linker in allowed_linkers:
        if any(mol.HasSubstructMatch(linker, atomIds=list(linker_atoms))):
            break
    else:
        return False, "Linker group between monomers not consistent with biflavonoids"
    
    # Check molecular weight and ring count
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for biflavonoid"
    
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5:
        return False, "Insufficient number of rings for biflavonoid"
    
    return True, "Contains at least two flavonoid monomers joined by a linker group consistent with biflavonoids"