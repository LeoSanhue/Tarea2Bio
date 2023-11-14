import requests
from bs4 import BeautifulSoup
import os
from collections import Counter
from Bio.PDB import PDBList, PDBParser
import py3Dmol
import tkinter as tk
from PIL import Image, ImageTk
import base64
import io

pdb_list = PDBList()
parser = PDBParser()


def descargar_proteina_pdb(codigo_pdb, archivo_html):
    url = f"https://www.rcsb.org/structure/{codigo_pdb}"
    response = requests.get(url)

    if response.status_code == 200:
        with open(archivo_html, "w", encoding="utf-8") as f:
            f.write(response.text)

        # print(f"\nInformación de {codigo_pdb} descargada con éxito en formato HTML.")
    else:
        print(f"Error al descargar información de {codigo_pdb}.")


def descargar_pdb_format(codigo_pdb, archivo_pdb):
    url = f"https://files.rcsb.org/download/{codigo_pdb.lower()}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        with open(archivo_pdb, "w", encoding="utf-8") as f:
            f.write(response.text)
        # print(f"\nArchivo PDB Format de {codigo_pdb} descargado con éxito.")
    else:
        print(f"Error al descargar el archivo PDB Format de {codigo_pdb}.")


def extraer_informacion(html_content, codigo_pdb, pdb_content, contenido_pdb):
    soup = BeautifulSoup(html_content, "html.parser")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(codigo_pdb, pdb_content)

    # 1) Nombre de la proteína y su código
    title = soup.title.text.strip()
    print(f"######################################################")
    print(f"\nNombre de la proteína y su código: {title}")

    # 2) Clasificaciones
    header_classification = soup.find("li", {"id": "header_classification"})
    if header_classification:
        classification_text = header_classification.strong.text.strip()
        print(f"{classification_text}")

    # 7) Información específica de la sección macromolecule-entityId-1-rowDescription
    macromolecule_info = soup.find(
        "tr", {"id": "macromolecule-entityId-1-rowDescription"}
    )
    if macromolecule_info:
        cells = macromolecule_info.find_all("td")

        nombre_proteina = cells[0].text.strip()
        cadenas = cells[1].find_all("a")
        cadena_info = ", ".join(
            [f"{cadena.text} ({cadena['href']})" for cadena in cadenas]
        )
        tamano_proteina = cells[2].text.strip()

        print(f"Nombre de la proteína: {nombre_proteina}")
        print(f"Cadenas: {cadena_info}")
        print(f"Tamaño de la proteína: {tamano_proteina}")
        contar_cadenas_repetidas(contenido_pdb, structure)


def contar_cadenas_repetidas(contenido_pdb, structure):
    # Dividir el contenido en líneas
    lineas = contenido_pdb.split("\n")

    # Inicializar un diccionario para contar las cadenas por tipo (A, B, C)
    cadenas_por_tipo = {"A": "", "B": "", "C": "", "D": ""}

    # Iterar sobre las líneas y actualizar las cadenas por tipo
    tipo_actual = None
    for linea in lineas:
        if linea.startswith("SEQRES"):
            tipo_actual = linea[11]  # Obtener el tipo de cadena (A, B, C)
            # Concatenar las cadenas eliminando "SEQRES" y espacios
            cadenas_por_tipo[tipo_actual] += linea[19:].replace(" ", "")

    # Contar frecuencia de cada cadena de 3 caracteres por tipo
    # tamaño = 0
    total_aa = 0
    for model in structure:
        for chain in model:
            amac = [residue.get_resname() for residue in chain if residue.id[0] == " "]
            largo = len(amac)
            total_aa += largo
            frecuencia = {aa: amac.count(aa) for aa in set(amac)}
            print(f"Tamaño Cadena {chain.id}: {largo}")
            print(f"Aminoacidos de la cadena {chain.id}: {frecuencia}\n")

    print(f"Total AA: {total_aa}")


def mostrar_estructura(contenido_pdb, nombre_proteina):
    viewer = py3Dmol.view()
    viewer.addModel(contenido_pdb, "pdb")
    viewer.setStyle({"cartoon": {"color": "spectrum"}})
    viewer.zoomTo()
    viewer.show()

    root = tk.Tk()
    root.title(nombre_proteina)

    # Crear una etiqueta para mostrar el nombre de la proteína
    label = tk.Label(root, text=f"Proteína: {nombre_proteina}")
    label.pack()

    # Crear un marco para contener la visualización 3D
    frame = tk.Frame(root)
    frame.pack()

    # Obtener el widget del visor 3D y embeberlo en el marco
    viewer_id = viewer.get_id()
    viewer_widget = tk.Frame(frame, width=800, height=400)
    viewer_widget.pack()
    viewer.tk.call(viewer_id, "embed", viewer_widget.winfo_id())

    root.mainloop()


# Elegir 5 proteínas de PDB
proteinas_pdb = ["6COX", "4INS", "1A0F", "1YMB", "1RHO"]

# Directorio para almacenar archivos
directorio_datos = "datos_proteinas_pdb"
if not os.path.exists(directorio_datos):
    os.makedirs(directorio_datos)

# Descargar información de cada proteína y extraer información
for proteina_pdb in proteinas_pdb:
    archivo_html = os.path.join(directorio_datos, f"{proteina_pdb}_info.html")
    pdb_file_path = pdb_list.retrieve_pdb_file(
        proteina_pdb, pdir=directorio_datos, file_format="pdb", overwrite=True
    )
    descargar_proteina_pdb(proteina_pdb, archivo_html)

    # Leer el contenido del archivo HTML
    with open(archivo_html, "r", encoding="utf-8") as f:
        contenido_html = f.read()

    # Llamar a la función para extraer información

    # Descargar el archivo PDB Format y contar cadenas repetidas
    archivo_pdb = os.path.join(directorio_datos, f"{proteina_pdb}.pdb")
    descargar_pdb_format(proteina_pdb, archivo_pdb)
    contenido_pdb = open(archivo_pdb, "r").read()
    extraer_informacion(contenido_html, proteina_pdb, archivo_pdb, contenido_pdb)
    # mostrar_estructura(contenido_pdb, proteina_pdb)
