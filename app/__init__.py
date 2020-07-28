from CGRdb import Molecule, load_schema
from CGRdb.database import MoleculeReaction
from CGRtools import MRVRead
from flask import Flask, render_template, request, redirect, url_for, abort
from io import BytesIO
from os import getenv
from json import loads
from pony.orm import db_session
from traceback import format_exc
from .watch import substructure, similarity, paths_of_synthesis_for_substructure_search, \
    paths_of_synthesis_for_target_molecule

app = Flask(__name__)

from .utils.list_converter import ListConverter

app.url_map.converters['list'] = ListConverter


@app.route("/")
@app.route("/check")
def check():
    return render_template('check.html')


@app.route("/index", methods=['POST'])
def index():
    try:
        forma = request.form
        if forma:
            if 'error' in forma:
                return render_template('index.html')
            elif 'pass' in forma:
                data = loads(request.form['pass'])
                if data['pass'] == 'password':
                    load_schema(getenv('DB_SCHEMA', 'name'), user=getenv('DB_USER', 'user'),
                                password=getenv('DB_PASS', 'password'), host=getenv('DB_HOST', 'host'),
                                database=getenv('DB_NAME', 'name'), port=getenv('DB_PORT', 'port'))
                    return render_template('index.html')
                else:
                    return render_template('check.html')
            else:
                return render_template('check.html')
    except Exception as e:
        return render_template('error.html', error=e), 402


@app.route("/post", methods=['POST'])
def post():
    try:
        if request.form:
            data = loads(request.form['data'])
            with BytesIO(data['source'].encode()) as f, MRVRead(f) as m:
                mol = next(m)
            mol.canonicalize()
            stages = int(data['stage'])
            available = bool(data['available'])
            with db_session:
                if data['substructure']:
                    return render_template("middle.html", molecules=substructure(mol), available=available,
                                           tree='hard_tree', stages=stages)
                else:
                    molecule = Molecule.find_structure(mol)
                    if molecule and MoleculeReaction.exists(molecule=molecule, is_product=True):
                        return redirect(
                            url_for('easy_tree', target=molecule.id, available=available, stages=stages))
                    else:
                        return render_template("middle.html", molecules=similarity(mol), available=available,
                                               tree='easy_tree', stages=stages)
    except Exception as e:
        print(format_exc())
        return render_template('error.html', error=e), 402


@app.route("/easy_tree/<int(min=1):target>/<available>/<int(min=1):stages>")
@db_session
def easy_tree(target, stages, available):
    try:
        target = Molecule[target]
    except:
        abort(404)
    try:
        if available == 'True':
            available = True
        else:
            available = False
        data = paths_of_synthesis_for_target_molecule(target, available, stages)
        return render_template('response.html', title='Molecule synthesis paths after searching molecule in database',
                               data=data)
    except Exception as e:
        print(format_exc())
        return render_template('error.html', error=e), 402


@app.route("/hard_tree/<int(min=1):target>/<available>/<list:reactions>/<int(min=1):stages>")
@db_session
def hard_tree(target, available, reactions, stages):
    try:
        target = Molecule[target]
    except:
        abort(404)
    if available == 'True':
        available = True
    else:
        available = False
    data = paths_of_synthesis_for_substructure_search(target, available, reactions, stages - 1)
    return render_template('response.html', title='Molecule synthesis paths after searching by substructure in database',
                           data=data)
