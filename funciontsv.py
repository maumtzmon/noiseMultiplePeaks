from pathlib import Path
import pandas as pd
import ast
import argparse


def resolve(df: pd.DataFrame, name: str) -> str:
    #Resolve column names in a DataFrame, ignoring case and spacing differences.
    key = name.strip().lower().replace(" ", "")
    for c in df.columns:
        if c.strip().lower().replace(" ", "") == key:
            return c
    raise KeyError(f"Column '{name}' not found. Got: {list(df.columns)}")

#Dado un RUNID valido y un MCMID valido, devuelve las listas de Data OK y Gain de ese row
#Usa "monitoring_DB.tsv" como default
#Si no encuentra el RUNID o MCMID, devuelve listas vacias
def listsFromTSV(path: str | None = None, run_ids=None, mcm_ids=None):

    if path is None:
        path = Path(__file__).with_name("monitoring_DB.tsv")
    else:
        path = Path(path)

    # Read TSV and normalize header spacing
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df = df.rename(columns={c: c.strip() for c in df.columns})

    run_col = resolve(df, "RUNID")
    mcm_col = resolve(df, "MCMID")
    dok_col = resolve(df, "Data OK")
    gain_col = resolve(df, "Gain") 

    # Filtra por RUNID 
    if run_ids is not None:
        if isinstance(run_ids, (int, float, str)):
            run_ids = [run_ids]
        try:
            run_ids = {int(str(r).strip()) for r in run_ids}
        except Exception:
            run_ids = {r for r in run_ids}
        run_series = pd.to_numeric(df[run_col].str.strip(), errors="coerce")
        df = df.loc[run_series.isin(run_ids)]

    # Filtra por MCMID
    if mcm_ids is not None:
        if isinstance(mcm_ids, (int, float, str)):
            mcm_ids = [mcm_ids]
        norm_mcm = {str(m).strip() for m in mcm_ids}
        df = df.loc[df[mcm_col].astype(str).str.strip().isin(norm_mcm)]

    # Build outputs
    def _parse_mcm(x: str):
        x = str(x).strip()
        try:
            return int(float(x))
        except Exception:
            return x

    listaMCMid = [_parse_mcm(x) for x in df[mcm_col].tolist()]

    #mete datos a la lista de Data OK
    listaDataOK: list[list[int]] = []
    for cell in df[dok_col].tolist():
        vals: list[int] = []
        for x in (ast.literal_eval(cell) if cell.startswith("[") and cell.endswith("]") else cell.strip("[]").split(",")):
            vals.append(int(float(x)))
        listaDataOK.append(vals)

    #mete datos a la lista de Gain
    listaGain: list[list[float]] = []
    for cell in df[gain_col].tolist():
        vals: list[float] = []
        for x in (ast.literal_eval(cell) if cell.startswith("[") and cell.endswith("]") else cell.strip("[]").split(",")):
                vals.append(float(str(x).replace(",", ".")))
        listaGain.append(vals)

    listaDataOK = [item for sublist in listaDataOK for item in sublist]
    listaGain = [item for sublist in listaGain for item in sublist]

    #copia el primer MCMID para que las listas queden del mismo tamaño
    if listaMCMid:
        listaMCMid = [listaMCMid[0]] * len(listaDataOK)
    else:
        listaMCMid = []

    return listaMCMid, listaDataOK, listaGain
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Lee el TSV de monitoreo y devuelve listas")
    parser.add_argument(
        "path",
        nargs="?",
        help="el path al TSV (default: monitoring_DB.tsv)",
    )
    parser.add_argument(
        "--run-ids",
        nargs="*",
        type=int,
    )

    parser.add_argument(
        "--mcm-ids",
        nargs="*",
    )
    args = parser.parse_args()

    mcmid, data_ok, gain = listsFromTSV(args.path, args.run_ids, args.mcm_ids)
    
    print(f"Path = {args.path or 'monitoring_DB.tsv'}")
    print(f"Run IDs = {args.run_ids}")
    print(f"MCM = {mcmid}")
    print(f"Tamaño MCM = {len(mcmid)}")
    print(f"Data OK = {data_ok}")
    print(f"Tamaño Data OK = {len(data_ok)}")
    print(f"Gain = {gain}")

