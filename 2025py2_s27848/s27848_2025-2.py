from Bio import Entrez, SeqIO
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import time
class NCBIRetriever:
    def __init__(self, email, api_key):
        self.email = email
        self.api_key = api_key
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'
        self.records_data = []
    def search_taxid(self, taxid, min_len, max_len):
        print(f"Szukam rekordów dla taxID: {taxid} z długością {min_len}-{max_len}")
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            tax_record = Entrez.read(handle)
            organism = tax_record[0]["ScientificName"]
            print(f"Organizm: {organism}")
            term = f"txid{taxid}[Organism] AND {min_len}:{max_len}[SLEN]"
            handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y")
            result = Entrez.read(handle)
            count = int(result["Count"])
            print(f"Znaleziono {count} rekordów.")
            if count == 0:
                return None
            self.webenv = result["WebEnv"]
            self.query_key = result["QueryKey"]
            self.count = count
            return count
        except Exception as e:
            print(f"Błąd podczas wyszukiwania: {e}")
            return None
    def fetch_records(self, max_to_fetch=None):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("Brak wyników wyszukiwania. Najpierw uruchom search_taxid().")
            return []
        max_records = min(self.count, max_to_fetch) if max_to_fetch else self.count
        print(f"Pobieranie {max_records} rekordów...")

        all_records = []
        for start in range(0, max_records, 500):
            print(f"Pobieranie: {start + 1} - {min(start + 500, max_records)}")
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=500,
                    webenv=self.webenv,
                    query_key=self.query_key
                )
                batch = list(SeqIO.parse(handle, "gb"))
                for record in batch:
                    all_records.append({
                        "Accession": record.id,
                        "Length": len(record.seq),
                        "Description": record.description
                    })
                time.sleep(1)
            except Exception as e:
                print(f"Błąd pobierania: {e}")
                continue
        self.records_data = all_records
        return all_records
    def save_csv(self, output_file):
        df = pd.DataFrame(self.records_data)
        df.to_csv(output_file, index=False)
        print(f"Zapisano dane do {output_file}")
    def plot_lengths(self, output_image):
        if not self.records_data:
            print("Brak danych do wykresu.")
            return
        df = pd.DataFrame(self.records_data)
        df_sorted = df.sort_values(by="Length", ascending=False)
        plt.figure(figsize=(14,6))
        plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o')
        plt.xticks(ticks=range(0,len(df_sorted), 10),labels=df_sorted["Accession"].iloc[::10],rotation=90, fontsize=6)
        plt.xlabel("Numer akcesyjny GenBank")
        plt.ylabel("Długość sekwencji")
        plt.title("Długości sekwencji GenBank")
        plt.grid(True, linestyle='--', linewidth=0.5)
        plt.tight_layout()
        plt.savefig(output_image)
        print(f"Wykres zapisany do {output_image}")
def main():
    parser = argparse.ArgumentParser(description="NCBI GenBank Retriever z CSV i wykresem")
    parser.add_argument("--email", required=True)
    parser.add_argument("--api-key", required=True)
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--min-len", type=int, default=200)
    parser.add_argument("--max-len", type=int, default=10000)
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()
    retriever = NCBIRetriever(args.email, args.api_key)
    count = retriever.search_taxid(args.taxid, args.min_len, args.max_len)
    if not count:
        print("Brak rekordów do pobrania.")
        return
    retriever.fetch_records(max_to_fetch=args.limit)
    retriever.save_csv(f"taxid_{args.taxid}_raport.csv")
    retriever.plot_lengths(f"taxid_{args.taxid}_wykres.png")
if __name__ == "__main__":
    main()