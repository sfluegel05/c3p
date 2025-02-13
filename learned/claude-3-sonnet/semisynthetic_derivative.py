"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: CHEBI:48911 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score

# Load and preprocess data (assuming you have a dataset of SMILES and labels)
data = [
    ('CC(=O)OC1=CC=CC=C1C(=O)O', 0),  # Non-derivative example
    ('CCCCN1CCC(=CC2=C1C=C(C=C2)OC)C(=O)OC', 1),  # Derivative example (dihydroergotamine)
    # Add more examples here...
]

def preprocess_data(data):
    X, y = [], []
    for smi, label in data:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            X.append(mol)
            y.append(label)
    return X, y

# Compute molecular descriptors
def compute_descriptors(mols):
    descriptors = []
    for mol in mols:
        desc = [
            Descriptors.MolWt(mol),
            Descriptors.NumRotatableBonds(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.TPSA(mol),
            # Add more descriptors if needed
        ]
        descriptors.append(desc)
    return descriptors

# Train a random forest model
def train_model(X, y):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    X_train_desc = compute_descriptors(X_train)
    X_test_desc = compute_descriptors(X_test)

    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train_desc, y_train)

    y_pred = model.predict(X_test_desc)
    f1 = f1_score(y_test, y_pred)
    return model, f1

# Classify a new molecule
def is_semisynthetic_derivative(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    desc = compute_descriptors([mol])
    prediction = model.predict(desc)[0]
    return bool(prediction), f"Classified as {'semisynthetic derivative' if prediction else 'not a semisynthetic derivative'}"

# Train the model
X, y = preprocess_data(data)
model, f1_score = train_model(X, y)
print(f"Model F1 score: {f1_score:.2f}")

# Example usage
smiles = "CCCCN1CCC(=CC2=C1C=C(C=C2)OC)C(=O)OC"  # Dihydroergotamine
prediction, reason = is_semisynthetic_derivative(smiles)
print(f"Prediction: {prediction}, Reason: {reason}")