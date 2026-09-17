###########截取fa文件头部和尾部2m区域
from Bio import SeqIO
import sys

def trim_fasta(input_file, output_header, output_tail, length=2000000):
    # 读取输入的fasta文件
    records = list(SeqIO.parse(input_file, 'fasta'))
    trimmed_records_header = []
    trimmed_records_tail = []

    # 遍历每个序列记录
    for record in records:
        # 提取头部和尾部的指定长度
        header_seq = record.seq[:length]
        tail_seq = record.seq[-length:]

        # 创建新的SeqRecord对象
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        header_record = SeqRecord(Seq(header_seq),
                                  id=record.id + "_header",
                                  description="Header 2M bases")
        tail_record = SeqRecord(Seq(tail_seq),
                                id=record.id + "_tail",
                                description="Tail 2M bases")

        # 添加到列表
        trimmed_records_header.append(header_record)
        trimmed_records_tail.append(tail_record)

    # 写入输出文件
    SeqIO.write(trimmed_records_header, output_header, "fasta")
    SeqIO.write(trimmed_records_tail, output_tail, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_fasta_file> <output_header_fasta> <output_tail_fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_header_fasta = sys.argv[2]
    output_tail_fasta = sys.argv[3]

    trim_fasta(input_fasta, output_header_fasta, output_tail_fasta)


