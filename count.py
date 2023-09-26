import os

# 获取当前目录
current_directory = os.getcwd()

# 初始化总文件大小为0
total_size = 0

# 遍历当前目录下的所有文件和文件夹
for root, dirs, files in os.walk(current_directory):
    for file in files:
        # 获取文件的完整路径
        file_path = os.path.join(root, file)
        
        # 排除以".root"结尾的文件
        if not (file.endswith(".root") or file.endswith(".dat") or "MuonHitDM" in root):
            # 获取文件大小并累加到总文件大小
            file_size = os.path.getsize(file_path)
            total_size += file_size

            # 判断文件大小是否超过10MB
            if file_size > 10 * 1024 * 1024:
                print("文件大小超过10MB的文件路径:", file_path)

# 打印总文件大小
print("总文件大小:", total_size, "字节")


# 转换为以MB或GB为单位的文件大小
if total_size < 1024 * 1024:
    # 小于1MB，输出以KB为单位的文件大小
    size_str = "{:.2f} KB".format(total_size / 1024)
elif total_size < 1024 * 1024 * 1024:
    # 小于1GB，输出以MB为单位的文件大小
    size_str = "{:.2f} MB".format(total_size / (1024 * 1024))
else:
    # 大于等于1GB，输出以GB为单位的文件大小
    size_str = "{:.2f} GB".format(total_size / (1024 * 1024 * 1024))

# 打印总文件大小
print("总文件大小:", size_str)
