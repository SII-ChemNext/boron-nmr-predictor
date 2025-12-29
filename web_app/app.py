"""Flask Web 应用 - 硼核磁预测系统"""

import os
import sys
from flask import Flask, render_template, request, jsonify, send_file, send_from_directory
from flask_cors import CORS
from datetime import datetime
import io
import csv

# 添加当前目录到 Python 路径
sys.path.insert(0, os.path.dirname(__file__))

from config import DevelopmentConfig
from core.predictor import BoronNMRPredictor
from database.models import init_db, save_prediction, get_history, get_prediction_by_id
from utils.validators import validate_input
from utils.exceptions import InvalidSMILESError, ValidationError, PredictionError

# 创建 Flask 应用
app = Flask(__name__, template_folder='templates', static_folder='static')
app.config.from_object(DevelopmentConfig)

# 启用 CORS
CORS(app)

# 全局预测器实例（应用启动时初始化）
predictor = None
_initialized = False


def initialize_app():
    """应用初始化"""
    global predictor, _initialized

    if _initialized:
        return

    _initialized = True

    print("\n" + "="*60)
    print("硼核磁预测系统 Web UI")
    print("="*60)

    # 初始化数据库
    try:
        init_db(app.config['DATABASE_PATH'])
    except Exception as e:
        print(f"✗ 数据库初始化失败: {e}")
        return

    # 创建图片输出目录
    os.makedirs(app.config['IMAGE_DIR'], exist_ok=True)

    # 加载模型
    print("\n正在加载模型...")
    try:
        predictor = BoronNMRPredictor(
            model_dir=app.config['MODEL_DIR'],
            device=app.config['DEVICE'],
            hidden_dim=app.config['HIDDEN_DIM'],
            dropout=app.config['DROPOUT']
        )
        print("✓ 模型加载完成\n")
    except Exception as e:
        print(f"✗ 模型加载失败: {e}\n")
        raise


# Flask 3.0+ 兼容的初始化方式
@app.before_request
def before_request():
    """在每个请求前执行"""
    initialize_app()


# ============================================================================
# 路由 - 页面
# ============================================================================

@app.route('/')
def index():
    """主页面"""
    solvents = app.config['SUPPORTED_SOLVENTS']
    return render_template('index.html', solvents=solvents)


@app.route('/history')
def history_page():
    """历史记录页面"""
    return render_template('history.html')


# ============================================================================
# Ketcher 编辑器路由
# ============================================================================

@app.route('/ketcher')
def ketcher_standalone():
    """Ketcher Standalone 编辑器 - 返回修改后的 index.html"""
    ketcher_path = os.path.join(
        os.path.dirname(__file__),
        'static/lib/ketcher/ketcher-standalone'
    )

    # 读取原始 index.html
    index_path = os.path.join(ketcher_path, 'index.html')
    with open(index_path, 'r', encoding='utf-8') as f:
        html = f.read()

    # 在 <head> 中添加 <base> 标签来修复相对路径
    # 将 ./static/ 替换为 /ketcher/static/
    html = html.replace('href="./static/', 'href="/ketcher/static/')
    html = html.replace('src="./static/', 'src="/ketcher/static/')

    # 修复其他相对路径
    html = html.replace('href="./', 'href="/ketcher/')
    html = html.replace('src="./', 'src="/ketcher/')

    return html


@app.route('/ketcher/<path:filename>')
def ketcher_static(filename):
    """Ketcher Standalone 的静态文件和资源"""
    ketcher_path = os.path.join(
        os.path.dirname(__file__),
        'static/lib/ketcher/ketcher-standalone'
    )
    return send_from_directory(ketcher_path, filename)


# ============================================================================
# API 路由 - 预测
# ============================================================================

@app.route('/api/predict', methods=['POST'])
def api_predict():
    """
    预测 API 端点

    请求 JSON:
    {
        "molecule_smiles": "OB(O)c1ccccc1",
        "solvent": "CDCl3"
    }

    响应 JSON:
    {
        "success": true,
        "prediction_id": 1,
        "canonical_smiles": "OB(O)c1ccccc1",
        "predictions": [
            {"atom_index": 1, "element": "B", "ppm": 28.53}
        ],
        "image_url": "/api/image/prediction_abc123.png",
        "num_borons": 1,
        "timestamp": "2025-12-24T10:30:00"
    }
    """
    try:
        if not predictor:
            return jsonify({
                'success': False,
                'error': '模型未加载，请稍后重试'
            }), 500

        # 1. 解析请求数据
        data = request.get_json()
        if not data:
            return jsonify({
                'success': False,
                'error': '请提供 JSON 数据'
            }), 400

        mol_smiles = data.get('molecule_smiles', '').strip()
        solvent_name = data.get('solvent', '').strip()

        # 2. 验证输入
        if not mol_smiles:
            return jsonify({
                'success': False,
                'error': 'molecule_smiles 不能为空'
            }), 400

        if not solvent_name:
            return jsonify({
                'success': False,
                'error': 'solvent 不能为空'
            }), 400

        try:
            canonical_smiles, solvent_smiles = validate_input(
                mol_smiles, solvent_name, app.config['SUPPORTED_SOLVENTS']
            )
        except (InvalidSMILESError, ValidationError) as e:
            return jsonify({
                'success': False,
                'error': str(e)
            }), 400

        # 3. 调用预测器
        try:
            result = predictor.predict(canonical_smiles, solvent_smiles)
        except PredictionError as e:
            return jsonify({
                'success': False,
                'error': f'预测失败: {str(e)}'
            }), 500

        # 4. 生成分子图片
        try:
            image_path = predictor.generate_molecule_image(
                result['mol_object'],
                [p['atom_index'] for p in result['predictions']],
                result['predictions'],
                app.config['IMAGE_DIR']
            )
            image_filename = os.path.basename(image_path)
        except Exception as e:
            app.logger.warning(f"分子图生成失败: {e}")
            image_filename = None

        # 5. 保存到数据库
        try:
            prediction_id = save_prediction(
                db_path=app.config['DATABASE_PATH'],
                mol_smiles=canonical_smiles,
                solvent_name=solvent_name,
                predictions=result['predictions'],
                image_path=image_filename if image_filename else ''
            )
        except Exception as e:
            app.logger.warning(f"数据库保存失败: {e}")
            prediction_id = None

        # 6. 返回结果
        return jsonify({
            'success': True,
            'prediction_id': prediction_id,
            'canonical_smiles': canonical_smiles,
            'predictions': result['predictions'],
            'image_url': f"/api/image/{image_filename}" if image_filename else None,
            'num_borons': result['num_borons'],
            'timestamp': datetime.now().isoformat()
        })

    except Exception as e:
        app.logger.error(f"预测过程异常: {e}", exc_info=True)
        return jsonify({
            'success': False,
            'error': f'服务器内部错误: {str(e)}'
        }), 500


# ============================================================================
# API 路由 - 历史记录
# ============================================================================

@app.route('/api/history')
def api_history():
    """
    获取历史记录 API

    查询参数:
    - limit: 返回的最大记录数（默认 50）

    响应 JSON:
    {
        "success": true,
        "records": [
            {
                "id": 1,
                "timestamp": "2025-12-24T10:30:00",
                "mol_smiles": "OB(O)c1ccccc1",
                "solvent_name": "CDCl3",
                "predictions": [{"atom_index": 1, "element": "B", "ppm": 28.53}],
                "image_path": "prediction_abc123.png"
            }
        ]
    }
    """
    try:
        limit = request.args.get('limit', app.config['HISTORY_LIMIT'], type=int)
        limit = min(limit, 200)  # 最多返回 200 条

        records = get_history(app.config['DATABASE_PATH'], limit=limit)

        return jsonify({
            'success': True,
            'records': records
        })

    except Exception as e:
        app.logger.error(f"获取历史记录失败: {e}")
        return jsonify({
            'success': False,
            'error': f'获取历史记录失败: {str(e)}'
        }), 500


@app.route('/api/prediction/<int:prediction_id>')
def api_get_prediction(prediction_id):
    """获取单个预测记录"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)
        return jsonify({
            'success': True,
            'record': record
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 404


# ============================================================================
# API 路由 - 文件服务
# ============================================================================

@app.route('/api/image/<filename>')
def api_serve_image(filename):
    """提供分子结构图"""
    try:
        # 安全检查：防止路径遍历
        if '..' in filename or '/' in filename:
            return jsonify({'error': '无效的文件名'}), 400

        img_path = os.path.join(app.config['IMAGE_DIR'], filename)

        if not os.path.exists(img_path):
            return jsonify({'error': '图片不存在'}), 404

        return send_file(img_path, mimetype='image/png')

    except Exception as e:
        app.logger.error(f"提供图片失败: {e}")
        return jsonify({'error': '文件服务错误'}), 500


# ============================================================================
# API 路由 - 下载
# ============================================================================

@app.route('/api/download/csv/<int:prediction_id>')
def api_download_csv(prediction_id):
    """下载 CSV 格式的预测结果"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)

        # 生成 CSV
        output = io.StringIO()
        writer = csv.writer(output)

        # 写入标题
        writer.writerow(['原子索引', '元素', '化学位移 (ppm)'])

        # 写入数据
        for pred in record['predictions']:
            writer.writerow([
                pred['atom_index'],
                pred['element'],
                f"{pred['ppm']:.2f}"
            ])

        # 写入元数据
        writer.writerow([])
        writer.writerow(['标准 SMILES', record['mol_smiles']])
        writer.writerow(['溶剂', record['solvent_name']])
        writer.writerow(['时间戳', record['timestamp']])

        # 返回文件
        output.seek(0)
        return send_file(
            io.BytesIO(output.getvalue().encode('utf-8')),
            mimetype='text/csv',
            as_attachment=True,
            download_name=f"prediction_{prediction_id}.csv"
        )

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/download/json/<int:prediction_id>')
def api_download_json(prediction_id):
    """下载 JSON 格式的预测结果"""
    try:
        record = get_prediction_by_id(app.config['DATABASE_PATH'], prediction_id)

        return send_file(
            io.BytesIO(
                jsonify(record).get_data(as_text=False).encode('utf-8')
            ),
            mimetype='application/json',
            as_attachment=True,
            download_name=f"prediction_{prediction_id}.json"
        )

    except Exception as e:
        return jsonify({'error': str(e)}), 500


# ============================================================================
# 错误处理
# ============================================================================

@app.errorhandler(404)
def not_found(error):
    """404 错误处理"""
    return jsonify({'error': '页面不存在'}), 404


@app.errorhandler(500)
def internal_error(error):
    """500 错误处理"""
    app.logger.error(f"内部错误: {error}")
    return jsonify({'error': '服务器内部错误'}), 500


@app.errorhandler(405)
def method_not_allowed(error):
    """405 方法不允许错误处理"""
    return jsonify({'error': '方法不允许'}), 405


# ============================================================================
# 应用启动
# ============================================================================

if __name__ == '__main__':
    # 确保在应用上下文中运行初始化
    with app.app_context():
        initialize_app()

    # 启动 Flask 应用
    app.run(
        host='0.0.0.0',
        port=5000,
        debug=app.config['DEBUG'],
        use_reloader=True
    )
