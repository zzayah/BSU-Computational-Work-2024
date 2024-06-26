from pvskiassess import Generation, Extraction

# Example Usage:
cluster_name = "your_cluster_name"
cluster_command = "your_cluster_command"
material = "SrCoO3_A.xyz"

generator = Generation(cluster_name, cluster_command)
generator.generate_files(material)
