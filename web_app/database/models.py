"""数据库模型 - 使用 SQLite"""

import sqlite3
import json
from datetime import datetime
import os


def init_db(db_path):
    """
    初始化数据库，创建必要的表

    Args:
        db_path (str): 数据库文件路径
    """
    # 创建目录如果不存在
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # 创建预测结果表
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
    print(f"✓ 数据库初始化完成: {db_path}")


def save_prediction(db_path, mol_smiles, solvent_name, predictions, image_path='', notes=''):
    """
    保存预测记录到数据库

    Args:
        db_path (str): 数据库文件路径
        mol_smiles (str): 分子 SMILES
        solvent_name (str): 溶剂名称
        predictions (list): 预测结果列表 [{'atom_index': int, 'element': str, 'ppm': float}, ...]
        image_path (str): 分子图片路径
        notes (str): 用户备注

    Returns:
        int: 插入的记录 ID

    Raises:
        Exception: 数据库操作失败
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # 转换预测结果为 JSON 字符串
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
        raise Exception(f"保存预测记录失败: {e}")


def get_history(db_path, limit=50):
    """
    从数据库获取历史记录

    Args:
        db_path (str): 数据库文件路径
        limit (int): 返回的最大记录数

    Returns:
        list: 预测记录列表

    Raises:
        Exception: 数据库操作失败
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
        raise Exception(f"获取历史记录失败: {e}")


def get_prediction_by_id(db_path, prediction_id):
    """
    根据 ID 获取单个预测记录

    Args:
        db_path (str): 数据库文件路径
        prediction_id (int): 预测记录 ID

    Returns:
        dict: 预测记录

    Raises:
        Exception: 记录不存在或数据库操作失败
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
            raise Exception(f"预测记录不存在: ID {prediction_id}")

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
        raise Exception(f"获取预测记录失败: {e}")


def delete_prediction(db_path, prediction_id):
    """
    删除预测记录

    Args:
        db_path (str): 数据库文件路径
        prediction_id (int): 预测记录 ID

    Returns:
        bool: 删除成功返回 True

    Raises:
        Exception: 数据库操作失败
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        cursor.execute('DELETE FROM predictions WHERE id = ?', (prediction_id,))
        conn.commit()
        conn.close()

        return True

    except Exception as e:
        raise Exception(f"删除预测记录失败: {e}")


def clear_old_predictions(db_path, days=30):
    """
    清理过期的预测记录

    Args:
        db_path (str): 数据库文件路径
        days (int): 保留多少天的记录

    Returns:
        int: 删除的记录数

    Raises:
        Exception: 数据库操作失败
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
        raise Exception(f"清理过期记录失败: {e}")
