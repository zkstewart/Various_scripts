import pyperclip

while True:
        # Obtain Apollo GFF3 annotation
        print("GO (Apollo annotation)")
        button=input()
        text=pyperclip.paste()
        lines=text.replace("\r", "").rstrip("\n").split("\n")
        rows =[]
        for line in lines:
                if line.startswith("##"):
                        continue
                rows.append(line.split("\t"))
        # Manually specify gene ID
        print("Enter ID for Apollo gene")
        geneID=input()
        out=[]
        exonCount=0
        cdsCount=0
        artifacts = []
        for row in rows:
                if row[2] == "gene":
                        details = "ID={0};Name=apollo_{0}".format(geneID)
                elif row[2] == "mRNA":
                        details="ID={0}.mrna1;Name=apollo_{0}.mrna1;Parent={0}".format(geneID)
                        # Extract mrnaID
                        mrnaID = [bit[3:] for bit in row[8].split(";") if bit.startswith("ID=")][0]
                elif row[2] == "exon":
                        exonCount += 1
                        details="Parent={0}.mrna1;ID={0}.mrna1.exon{1};Name=apollo_{0}.mrna1".format(geneID, exonCount)
                elif row[2] == "CDS":
                        cdsCount += 1
                        details="Parent={0}.mrna1;ID={0}.mrna1.cds{1};Name=apollo_{0}.mrna1".format(geneID, cdsCount)
                elif "artifact" in row[2]:
                        artifactType = row[2].split("_")[0]
                        if "residues" in row[8]:
                                residues = ";residues={0}".format([bit.split("=")[1] for bit in row[8].split(";") if bit.startswith("residues=")][0])
                        else:
                                residues = ""
                        details="ID={0}_for_{1}.mrna1;Name=apollo_{0}_for_{1}.mrna1{2}".format(artifactType, geneID, residues)
                row[8] = details
                row[1] = "apollo"
                if "artifact" not in row[2]:
                        out.append("\t".join(row))
                else:
                        artifacts.append("\t".join(row))
        # Obtain protein
        print("Enter the protein sequence for gene")
        protein=input()
        if protein == "":
                protein = pyperclip.paste()
                protein = protein.rstrip("\r\n ")
                protein = protein.replace("\r", "").replace("\n", "")
        # Format output
        out.insert(0, "# APOLLO ANNOTATION: {0}.mrna1 manual model build".format(geneID))
        out.append("#PROT {0}.mrna1 {0}\t{1}".format(geneID, protein))
        out += artifacts
        pyperclip.copy("\n".join(out))
