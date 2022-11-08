const { ipcRenderer } = require("electron");
const math = require("mathjs");

function run() {
var n = document.getElementById('n').value;
var sum_Xi = document.getElementById('sum_Xi').value;
var sum_Zi = document.getElementById('sum_Zi').value;
var sum_Yi = document.getElementById('sum_Yi').value;
var sum_Xi_sqr = document.getElementById('sum_Xi_sqr').value;
var sum_Zi_sqr = document.getElementById('sum_Zi_sqr').value;
var sum_Yi_sqr = document.getElementById('sum_Yi_sqr').value;
var sum_XiZi = document.getElementById('sum_XiZi').value;
var sum_YiZi = document.getElementById('sum_YiZi').value;
var sum_XiYi = document.getElementById('sum_XiYi').value;
var k = document.getElementById('k').value;

const matrix1 = [[n, sum_Xi ,sum_Zi ], [sum_Xi ,sum_Xi_sqr ,sum_XiZi ], [sum_Zi ,sum_XiZi , sum_Zi_sqr ]]
const matrix_XtX = math.matrix(matrix1)    

const matrix2 = [sum_Yi , sum_XiYi , sum_YiZi ]
const matrix_XtY = math.matrix(matrix2)

const det = math.det(matrix_XtX)

const inverse = (m) => {
    return math.inv(m);
  };
  const matrix_XtX_inverse = inverse(matrix_XtX);

  const betaMu = math.multiply(matrix_XtX_inverse, matrix_XtY)

  document.getElementById('result1').innerHTML = betaMu


var sum_resid_sqr = sum_Yi_sqr - (betaMu[0]*sum_Yi + betaMu[1]*sum_XiYi + betaMu[2]*sum_YiZi);
var sigmaMu_sqr = sum_resid_sqr/(n-k);
const khoảng_tin_cậy_betaMu = []
const standardError_array = []

for (j = 1; j < betaMu.length; j++){
    var variance_betaMu = sigmaMu_sqr*(matrix_XtX_inverse[j][j])
    var standardError = Math.sqrt(variance_betaMu)
    standardError_array.push(standardError)
    var KTC_left = betaMu - 2.626*standardError
    var KTC_right = betaMu + 2.626*standardError
    khoảng_tin_cậy_betaMu.push(KTC_left + '-' + KTC_right)
}
document.getElementById('result2').innerHTML = khoảng_tin_cậy_betaMu
}