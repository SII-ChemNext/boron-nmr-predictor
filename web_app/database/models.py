"""Database models - using SQLite"""

import sqlite3
import json
from datetime import datetime
import os


def init_db(db_path):
    """
    Initialize the database and create required tables.

    Args:
        db_path (str): path to the database file
    """
    # Create directory if it does not exist
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create predictions table
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS predictions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
            mol_smiles TEXT NOT NULL,
            solvent_name TEXT NOT NULL,
            predictions_json TEXT NOT NULL,
            image_path TEXT,
            notes TEXT
        )
    ''')

    conn.commit()
    conn.close()
    print(f"✓ Database initialized: {db_path}")


def save_prediction(db_path, mol_smiles, solvent_name, predictions, image_path='', notes=''):
    """
    Save a prediction record to the database.

    Args:
        db_path (str): path to the database file
        mol_smiles (str): molecule SMILES
        solvent_name (str): solvent name
        predictions (list): list of prediction results [{'atom_index': int, 'element': str, 'ppm': float}, ...]
        image_path (str): path to the molecule image
        notes (str): user notes

    Returns:
        int: ID of the inserted record

    Raises:
        Exception: if database operation fails
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Serialize prediction results to JSON string
        predictions_json = json.dumps(predictions, ensure_ascii=False)

        cursor.execute('''
            INSERT INTO predictions (mol_smiles, solvent_name, predictions_json, image_path, notes)
            VALUES (?, ?, ?, ?, ?)
        ''', (mol_smiles, solvent_name, predictions_json, image_path, notes))

        prediction_id = cursor.lastrowid
        conn.commit()
        conn.close()

        return prediction_id

    except Exception as e:
        raise Exception(f"Failed to save prediction record: {e}")


def get_history(db_path, limit=50):
    """
    Retrieve prediction history from the database.

    Args:
        db_path (str): path to the database file
        limit (int): maximum number of records to return

    Returns:
        list: list of prediction records

    Raises:
        Exception: if database operation fails
    """
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        cursor.execute('''
            SELECT * FROM predictions
            ORDER BY timestamp DESC
            LIMIT ?
        ''', (limit,))

        rows = cursor.fetchall()
        conn.close()

        records = []
        for row in rows:
            records.append({
                'id': row['id'],
                'timestamp': row['timestamp'],
                'mol_smiles': row['mol_smiles'],
                'solvent_name': row['solvent_name'],
                'predictions': json.loads(row['predictions_json']),
                'image_path': row['image_path'],
                'notes': row['notes']
            })

        return records

    except Exception as e:
        raise Exception(f"Failed to retrieve history: {e}")


def get_prediction_by_id(db_path, prediction_id):
    """
    Retrieve a single prediction record by ID.

    Args:
        db_path (str): path to the database file
        prediction_id (int): prediction record ID

    Returns:
        dict: prediction record

    Raises:
        Exception: if record does not exist or database operation fails
    """
    try:
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        cursor.execute('''
            SELECT * FROM predictions
            WHERE id = ?
        ''', (prediction_id,))

        row = cursor.fetchone()
        conn.close()

        if row is None:
            raise Exception(f"Prediction record not found: ID {prediction_id}")

        return {
            'id': row['id'],
            'timestamp': row['timestamp'],
            'mol_smiles': row['mol_smiles'],
            'solvent_name': row['solvent_name'],
            'predictions': json.loads(row['predictions_json']),
            'image_path': row['image_path'],
            'notes': row['notes']
        }

    except Exception as e:
        raise Exception(f"Failed to retrieve prediction record: {e}")


def delete_prediction(db_path, prediction_id):
    """
    Delete a prediction record.

    Args:
        db_path (str): path to the database file
        prediction_id (int): prediction record ID

    Returns:
        bool: True if deletion was successful

    Raises:
        Exception: if database operation fails
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        cursor.execute('DELETE FROM predictions WHERE id = ?', (prediction_id,))
        conn.commit()
        conn.close()

        return True

    except Exception as e:
        raise Exception(f"Failed to delete prediction record: {e}")


def clear_old_predictions(db_path, days=30):
    """
    Remove expired prediction records.

    Args:
        db_path (str): path to the database file
        days (int): number of days of records to retain

    Returns:
        int: number of deleted records

    Raises:
        Exception: if database operation fails
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        cursor.execute('''
            DELETE FROM predictions
            WHERE datetime(timestamp) < datetime('now', '-' || ? || ' days')
        ''', (days,))

        deleted_count = cursor.rowcount
        conn.commit()
        conn.close()

        return deleted_count

    except Exception as e:
        raise Exception(f"Failed to clear expired records: {e}")
