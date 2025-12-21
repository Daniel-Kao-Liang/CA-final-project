#!/bin/bash

# 检查是否传入了<student id>
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <student id>"
    exit 1
fi

STUDENT_ID=$1

# 创建<student id>文件夹
if [ -d "$STUDENT_ID" ]; then
    echo "Directory $STUDENT_ID already exists. Deleting it."
    rm -rf "$STUDENT_ID"
fi
mkdir "$STUDENT_ID"

# 复制当前目录的所有文件到<student id>文件夹中
for file in *; do
    if [ "$file" != "$0" ] && [ "$file" != "$STUDENT_ID" ]; then
        cp -r "$file" "$STUDENT_ID/"
    fi
done

# 压缩为<student id>.tar.gz
tar -czf "$STUDENT_ID.tar.gz" "$STUDENT_ID"

# 提示压缩完成
echo "Files have been compressed into $STUDENT_ID.tar.gz."
