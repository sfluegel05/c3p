"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:64211 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids include endocannabinoids (ethanolamides, glycerol esters) and
    phytocannabinoids (dibenzopyran derivatives), characterized by specific oxygen-containing groups.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check ethanolamide group (N-linked to carbonyl and ethanolamine)
    ethanolamide_pattern = Chem.MolFromSmarts("[CX3](=O)-[NH]-[CH2]-[CH2]-O")
    if mol.HasSubstructMatch(ethanolamide_pattern):
        return True, "Ethanolamide group present (endocannabinoid)"

    # Check glycerol ester/ether with long chain (mono- or di-substituted)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH](O[!H0])[CH2]")
    if mol.HasSubstructMatches(glycerol_pattern):
        # Look for at least one long chain (>=16 carbons) connected via ester/ether
        chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4]")
        for match in mol.GetSubstructMatches(glycerol_pattern):
            env_atoms = set()
            for atom_idx in match:
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, 6, atom_idx)
                env_atoms.update(env)
            chain_matches = len(mol.GetSubstructMatches(chain_pattern, uniquify=False, maxMatches=20))
            if chain_matches >= 4:  # Approximate long chain presence
                return True, "Glycerol derivative with long chain"

    # Check dibenzopyran core (phytocannabinoids like THC)
    dibenzopyran_core = Chem.MolFromSmarts("[c]1[c][c][c]([C@@H]2[C@H](CC[C@@H]2O)C)[c]([O])[c]1")
    if mol.HasSubstructMatch(dibenzopyran_core):
        return True, "Dibenzopyran core detected"

    # Check synthetic scaffolds (indole/indazole + carbonyl linker)
    indole_pattern = Chem.MolFromSmarts("[nH]1c([CX3]=O)ccc2ccccc12")  # Indole carbonyl
    indazole_pattern = Chem.MolFromSmarts("[nH]1c([CX3]=O)nnc2ccccc12")  # Indazole carbonyl
    if mol.HasSubstructMatch(indole_pattern) or mol.HasSubstructMatch(indazole_pattern):
        return True, "Synthetic cannabinoid scaffold detected"

    # Oxygen in heterocyclic ring (e.g., tetrahydrocannabinol variants)
    oxygen_rings = False
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        if len(ring) < 5:
            continue
        oxygen_present = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if oxygen_present and any(mol.GetAtomWithIdx(a).GetDegree() >= 2 for a in ring):
            oxygen_rings = True
            break
    if oxygen_rings:
        return True, "Oxygen in heterocyclic ring system"

    # Check for polyunsaturated long chain (>=18 carbons) with oxygen groups
    unsaturated_chain = Chem.MolFromSmarts("C/C=C/C/C=C/C/C=C/C/C=C/CCCC(=O)")
    if mol.HasSubstructMatch(unsaturated_chain):
        return True, "Polyunsaturated long chain with carbonyl group"

    return False, "No characteristic cannabinoid features detected"