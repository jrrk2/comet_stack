#ifndef IMAGE_READER_H
#define IMAGE_READER_H

#include <QString>
#include <QVector>
#include <QStringList>
#include <memory>

// Forward declarations
class ImageReaderPrivate;

struct ImageData {
    int width = 0;
    int height = 0;
    int channels = 0;
    QVector<float> pixels;
    QString colorSpace;
    QString format;
    QStringList metadata;
    
    bool isValid() const { 
        return width > 0 && height > 0 && channels > 0 && !pixels.isEmpty(); 
    }
    
    void clear() {
        width = height = channels = 0;
        pixels.clear();
        colorSpace.clear();
        format.clear();
        metadata.clear();
    }
};

class ImageReader
{
public:
    explicit ImageReader();
    ~ImageReader();

    // File reading
    bool readFile(const QString& filePath);
    
    // Data access
    const ImageData& imageData() const;
    bool hasImage() const;
    void clear();
    void setImageData(const ImageData& imageData);
    
    // Error handling
    QString lastError() const;
    
    // Supported formats
    static QStringList supportedFormats();
    static QString formatFilter();
    
    // File type detection
    static QString detectFileType(const QString& filePath);
    
private:
    std::unique_ptr<ImageReaderPrivate> d;
    
    // Non-copyable
    ImageReader(const ImageReader&) = delete;
    ImageReader& operator=(const ImageReader&) = delete;
};

#endif // IMAGE_READER_H
