/**
 *
 * @param {*} content :Blast+ outfmt 0
 * author:zpliu
 * Date:2020/04/08
 */

// var fs = require('fs')
// fs.readFile('./1', 'utf-8', function (err, data) {
//   console.log(data)
// })
// text = fs.readFileSync('./1', 'utf-8')
// parseBlastText(text)
// console.log(text)
function parseBlastText(content) {
  var lines = content.split('\n')
  var BlastTextJson = []
  var queryIndexArray = []
  /*
   *获取query 基因编号所在行
   *为了方便遍历，将最后一行的行号加入数组
   */
  for (var i = 0; i < lines.length; i++) {
    if (lines[i].startsWith('Query= ')) {
      queryIndexArray.push(i)
    }
  }
  queryIndexArray.push(lines.length - 1)
  // console.log(lines[queryIndexArray[2]])
  /*
   * 遍历每个query所在的区域
   */
  for (var index = 0; index < queryIndexArray.length - 1; index++) {
    var query = ''
    var queryLength = ''
    var alignments = []
    for (var i = queryIndexArray[index]; i < queryIndexArray[index + 1]; i++) {
      if (lines[i].startsWith('Query= ')) {
        query = lines[i].split('Query= ')[1]
        queryLength = lines[i + 2].split('Length=')[1]
      }
      /*
       * 获取subjec所在行编号
       */
      if (lines[i].startsWith('>')) {
        var line1 = lines[i].split('>')[1]
        var line2 = ''
        var currentLine = i
        while (line2 == '') {
          currentLine = currentLine + 1
          if (lines[currentLine].startsWith('Length=')) {
            line2 = lines[currentLine]
          } else {
            line1 += lines[currentLine]
          }
        }
        var subjectName = line1.split(/\s+/)[0]
        var subjectDescription = line1.split(subjectName)[1].replace(/\s+/, '')
        var subjectLength = line2.split('=')[1]
        var score = lines[currentLine + 2]
          .split(',')[0]
          .replace(/\s\s+/g, ' ')
          .split(' ')[3]
        var eValue = lines[currentLine + 2].split(',')[1].split(' ')[4]
        var identities = lines[currentLine + 3]
          .split(',')[0]
          .split('(')[1]
          .substr(
            0,
            lines[currentLine + 3].split(',')[0].split('(')[1].length - 2
          )
        if (lines[0].startsWith('BLASTN')) {
          var positives = 'N/A'
          var gaps = lines[currentLine + 3]
            .split(',')[1]
            .split('(')[1]
            .substr(
              0,
              lines[currentLine + 3].split(',')[1].split('(')[1].length - 2
            )
        } else {
          var positives = lines[currentLine + 3]
            .split(',')[1]
            .split('(')[1]
            .substr(
              0,
              lines[currentLine + 3].split(',')[1].split('(')[1].length - 2
            )
          var gaps = lines[currentLine + 3]
            .split(',')[2]
            .split('(')[1]
            .substr(
              0,
              lines[currentLine + 3].split(',')[2].split('(')[1].length - 2
            )
        }
        /*
         * 获取序列比对情况的范围
         */
        var pairSequenIndex = currentLine + 4
        while (lines[pairSequenIndex] + lines[pairSequenIndex + 1] != '') {
          pairSequenIndex += 1
        }
        var pairSequence = ''
        for (var i = currentLine + 5; i < pairSequenIndex; i = i + 4) {
          pairSequence +=
            lines[i] +
            '\n' +
            lines[i + 1] +
            '\n' +
            lines[i + 2] +
            '\n' +
            lines[i + 3] +
            '\n'
        }
        // var queryStart = ''
        // var queryEnd = ''
        // var subjectStart = ''
        // var subjectEnd = ''
        // var sequence = ''
        // for (var i = currentLine + 5; i < pairSequenIndex; i = i + 4) {
        //   queryStart = lines[i].split(/\s+/)[1]
        //   queryEnd = lines[i].split(/\s+/)[3]
        //   subjectStart = lines[i].split(/\s+/)[1]
        //   subjectEnd = lines[i].split(/\s+/)[3]
        //   sequence =
        //     lines[i].split(/\s+/)[2] +
        //     '<br />' +
        //     lines[i + 1].slice(12) +
        //     '<br/>' +
        //     lines[i + 2].split(/\s+/)[2]
        //   pairSequence.push({
        //     queryStart,
        //     queryEnd,
        //     subjectStart,
        //     subjectEnd,
        //     sequence
        //   })
        // }
        /*
         *将该query比对信息,封装成Array
         */
        alignments.push({
          subjectName,
          subjectLength,
          subjectDescription,
          score,
          eValue,
          positives,
          pairSequence,
          gaps,
          identities
        })
      }
    }
    BlastTextJson.push({
      query,
      queryLength,
      alignments
    })
  }
  return BlastTextJson
}
module.exports = parseBlastText
